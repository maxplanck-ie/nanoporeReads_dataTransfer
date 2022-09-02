import sys
import os
import glob
import pandas as pd
import numpy as np
from rich import print
from importlib.metadata import version
import yaml
import json
import subprocess as sp
from npr.communication import send_email
from npr.snakehelper import glob2reports
from npr.snakehelper import getfoot

def find_new_flowcell(config):
    """
    look for new flowcells that are flagged with `transfer.final`.
    """
    # set base path
    base_path = os.path.join(
        config["paths"]["baseDir"]
    )
    # glob dirs with a transfer flag.
    dirs = glob.glob(
        os.path.join(
            base_path, "*/transfer.final"
        )
    )
    # Iterate over dirs.
    for dir in dirs:
         # abs path to flowcell
        flowcell = os.path.dirname(dir)
        # trigger when no flowcell folder in output directory.
        if not os.path.exists(
            os.path.join(
                config["paths"]["outputDir"],
                os.path.basename(flowcell),
                'analysis.done'
            )
        ):
            msg = "New flowcell found: {} \n".format(
                os.path.basename(flowcell)
            )
            if not os.path.isfile(
                os.path.join(flowcell, 'SampleSheet.csv')
            ):
                msg += "No SampleSheet.csv file.. Exiting.\n"
                send_email(
                    msg,
                    version("npr"),
                    os.path.basename(flowcell),
                    config
                )
                sys.exit("no sampleSheet.")
            else:
                msg += "SampleSheet.csv file found. Start processing..\n"
                # send_email(
                #    msg,
                #    version("npr"),
                #    os.path.basename(flowcell),
                #    config
                # )
                config["input"] = {
                    'name': os.path.basename(flowcell)
                }
                return (os.path.basename(flowcell), msg)
    return (None, None)

def read_flowcell_info(config, info_dict):
    """
     - set paths
     - copy json, summaries & samplesheet
     - parse json for flowcell, kit & barcoding
     - infer model from flowcell + kit
     - convert fast5 to pod5
    """
    flowcell = config["input"]["name"]
    base_path = os.path.join(
        config["paths"]["baseDir"],
        flowcell
    )
    info_dict["base_path"] = base_path
    flowcell_path = os.path.join(
        config["paths"]["outputDir"],
        flowcell
    )
    if not os.path.exists(flowcell_path):
        os.mkdir(
            flowcell_path
        )
    info_dict["flowcell_path"] = flowcell_path
    # Copy over reports.
    for globStr in [
        '*json',
        '*html',
        "*txt",
        'SampleSheet.csv'
    ]:
        glob2reports(
            globStr,
            base_path,
            flowcell_path
        )
    json_file = glob.glob(
        os.path.join(
            flowcell_path,
            'reports',
            '*json'
        )
    )
    with open(json_file[0]) as f: # assume only 1 file
        jsondata = json.load(f)
    info_dict["flowcell"] = jsondata['protocol_run_info']['meta_info']['tags']['flow cell']['string_value']
    info_dict["kit"] = jsondata['protocol_run_info']['meta_info']['tags']['kit']['string_value']
    info_dict['barcoding'] = bool(jsondata['protocol_run_info']['meta_info']['tags']['barcoding']['bool_value'])
    print("flowcell = {}".format(info_dict["flowcell"]))
    print("kit = {}".format(info_dict["kit"]))
    modeldic = yaml.safe_load(
        open(config['guppy_basecaller']['models'])
    )
    info_dict['model'] = modeldic[info_dict['flowcell']][info_dict['kit']]
    print('model = {}'.format(info_dict['model']))
    poddir = os.path.join(
        flowcell_path,
        'pod5'
    )
    if not os.path.exists(poddir):
        os.mkdir(poddir)
    # podracing time
    if not os.path.exists(
        os.path.join(poddir, 'output.pod5')
    ):
        print("[bold green]Converting fast5-pod5[/bold green]")
        podracercmd = [
            'pod5-convert-from-fast5',
            base_path,
            poddir,
            '--recursive',
            '--force-overwrite',
            '--active-readers',
            '10'
        ]
        sp.check_output(podracercmd)
    else:
        print("[bold green]pod5 exists, not converting[/bold green]")
    info_dict['poddir'] = poddir
    # Check sizes for ratios.
    original_foot = getfoot(
        os.path.join(
            base_path,
            'fast5_pass'
        )
    )
    original_foot += getfoot(
        os.path.join(
            base_path,
            'fast5_fail'
        )
    )

    new_foot = getfoot(poddir)
    print(original_foot)
    print(new_foot)
    info_dict['pod5 compression'] = round(new_foot/original_foot, 2)
    return(info_dict)

def read_samplesheet(config):
    """
    read samplesheet
    """
    sample_sheet = pd.read_csv(
        os.path.join(
            config["info_dict"]['flowcell_path'],
            'reports',
            'SampleSheet.csv'
        ),
        sep = ",",
        skiprows=[0]
    )
    bc_kit = 'bc' if config['info_dict']['barcoding'] else 'no_bc'
    # parkour doesn't like 'no indices', there's thus an index list with 'No_index1, No_index2, ...'
    sample_sheet['I7_Index_ID'] = sample_sheet['I7_Index_ID'].str.replace('No_index*','no_bc', regex = True)
    if any(sample_sheet['I7_Index_ID'].str.contains('no_bc')):
       bc_kit = "no_bc"
    if bc_kit == 'bc':
       bc_kit = np.unique(sample_sheet["Description"].values)[0]
       start = bc_kit.find("(") + len("(")
       end = bc_kit.find(")")
       bc_kit = bc_kit[start:end]
    data=dict()
    data['projects'] = []
    data['samples'] = []
    for index, row in sample_sheet.iterrows():
        assert(row["Sample_ID"] not in data.keys())
        if row["Sample_Project"] not in data['projects']:
            data['projects'].append(row["Sample_Project"])
        if row['Sample_ID'] not in data['samples']:
            data['samples'].append(row['Sample_ID'])
        data[row["Sample_ID"]] = dict({"Sample_Name": row["Sample_Name"],
                                       "Sample_Project": row["Sample_Project"],
                                       "barcode_kits": bc_kit,
                                       "index_id": row["I7_Index_ID"].replace(
                                        'BP', 'barcode'
                                       ),
                                       "Sample_ID": row["Sample_ID"]})
    print("[green] Barcode kit determined as: {}".format(bc_kit))
    return bc_kit, data
