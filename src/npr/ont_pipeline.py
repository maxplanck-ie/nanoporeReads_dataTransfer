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
    if json_file:
        with open(json_file[0]) as f: # assume only 1 file
            jsondata = json.load(f)
        info_dict["flowcell"] = jsondata['protocol_run_info']['meta_info']['tags']['flow cell']['string_value']
        info_dict["kit"] = jsondata['protocol_run_info']['meta_info']['tags']['kit']['string_value']
        info_dict['barcoding'] = bool(jsondata['protocol_run_info']['meta_info']['tags']['barcoding']['bool_value'])
        # double check args. This needs a cleaner solution.
        for rg in jsondata['protocol_run_info']['args']:
            if rg == '--barcoding' and not info_dict['barcoding']:
                print("Json bool for barcoding wrong ! Override !")
                info_dict['barcoding'] = True
            if rg.startswith('barcoding_kits'):
                start = rg.find('[') + 1
                stop = rg.find(']')
                info_dict['barcode_kit'] = rg[start:stop].replace('\"', '')
        if 'barcode_kit' not in info_dict:
            info_dict['barcode_kit'] = 'no_bc'
    else:
        # try to get md file.
        finsum = glob.glob(
            os.path.join(
                base_path,
                'final_summary*txt'
            )
        )[0]
        if finsum:
            for line in open(finsum, 'r'):
                if line.strip().startswith('protocol='):
                    info_dict['flowcell'] = line.strip().split(':')[-2]
                    info_dict['kit'] = line.strip().split(':')[-1]
                    print("assuming no barcode.")
                    info_dict['barcode_kit'] = 'no_bc'
                    info_dict['barcoding'] = False
        else:
            sys.exit("no json file, no final summary txt file found. exiting.")
    print("flowcell = {}".format(info_dict["flowcell"]))
    print("kit = {}".format(info_dict["kit"]))
    modeldic = yaml.safe_load(
        open(config['guppy_basecaller']['models'])
    )
    info_dict['model'] = modeldic[info_dict['flowcell']][info_dict['kit']]
    print('model = {}'.format(info_dict['model']))
    poddirpass = os.path.join(
        flowcell_path,
        'pod5_pass'
    )
    poddirfail = os.path.join(
        flowcell_path,
        'pod5_fail'
    )
    info_dict['poddirpass'] = poddirpass
    info_dict['poddirfail'] = poddirfail
    return (info_dict)

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

    sample_sheet['I7_Index_ID'] = sample_sheet['I7_Index_ID'].str.replace('No_index*','no_bc', regex = True)


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
                                       "barcode_kits": config["info_dict"]['barcode_kit'],
                                       "index_id": row["I7_Index_ID"].replace(
                                        'BP', 'barcode'
                                       ).replace(
                                        'NB', 'barcode'
                                       ),
                                       "Sample_ID": row["Sample_ID"]})
    print("[green] Barcode kit determined as: {}".format(
        config["info_dict"]['barcode_kit']
    ))
    return config["info_dict"]['barcode_kit'], data
