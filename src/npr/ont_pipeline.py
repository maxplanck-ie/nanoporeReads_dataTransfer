import sys
import os
import glob
import shutil
import pandas as pd
import numpy as np
from rich import print
from importlib.metadata import version
from npr.communication import send_email
import yaml
import json
import subprocess as sp

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
    Check the flowcell path to find the info needed for base calling
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
    # json files
    jsonfiles = glob.glob(
        os.path.join(
            base_path,
            '*json'
        )
    )
    if jsonfiles == []:
        sys.exit("no json file found..")
    shutil.copy(
        jsonfiles[0],
        os.path.join(
            flowcell_path,
            os.path.basename(jsonfiles[0])
        )
    )
    with open(jsonfiles[0]) as f: #assume there is always only 1 json file..
        jsondata = json.load(f)
    info_dict["flowcell"] = jsondata['protocol_run_info']['meta_info']['tags']['flow cell']['string_value']
    info_dict["kit"] = jsondata['protocol_run_info']['meta_info']['tags']['kit']['string_value']
    info_dict['barcoding'] = bool(jsondata['protocol_run_info']['meta_info']['tags']['barcoding']['bool_value'])
    print("flowcell = {}".format(info_dict["flowcell"]))
    print("kit = {}".format(info_dict["kit"]))
    # summary files.
    summaryfiles = glob.glob(
        os.path.join(
            base_path,
            '*summary*'
        )
    )
    if summaryfiles == []:
        sys.exit('no summary files found..')
    for s in summaryfiles:
        shutil.copy(
            s,
            os.path.join(
                flowcell_path,
                os.path.basename(s)
            )
        )
    # samplesheet.
    if not os.path.exists(
        os.path.join(base_path, 'SampleSheet.csv')
    ):
        sys.exit("no samplesheet found.")
    else:
        shutil.copy(
            os.path.join(base_path, 'SampleSheet.csv'),
            os.path.join(flowcell_path, 'SampleSheet.csv')
        )
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
        podracercmd = [
            'pod5-convert-from-fast5',
            base_path,
            poddir,
            '--recursive',
            '--force-overwrite',
            '--active-readers',
            '20'
        ]
        sp.check_output(podracercmd)
    info_dict['poddir'] = poddir
    # Check sizes for ratios.
    original_foot = 0
    for f5d in glob.glob(
        os.path.join(
            base_path,
            'fast5*'
        )
    ):
        for f in os.listdir(f5d):
            original_foot += os.path.getsize(
                os.path.join(f5d, f)
            )
    new_foot = os.path.getsize(
        os.path.join(
            poddir,
            'output.pod5'
        )
    )
    info_dict['pod5 compression'] = round(new_foot/original_foot, 2)
    return(info_dict)


def read_samplesheet(config):
    """
    read samplesheet
    """
    sample_sheet = pd.read_csv(
        config["info_dict"]["flowcell_path"]+"/SampleSheet.csv",
        sep = ",",
        skiprows=[0]
    )
    bc_kit = 'bc' if config['info_dict']['barcoding'] else 'no_bc'
    # parkour doesn't like 'no indices', there's thus an index list with 'No_index1, No_index2, ...'
    #sample_sheet['I7_Index_ID'] = sample_sheet['I7_Index_ID'].str.replace('No_index*','no_bc', regex = True)
    #if any(sample_sheet['I7_Index_ID'].str.contains('no_bc')):
    #    bc_kit = "no_bc"
    if bc_kit == 'bc':
       bc_kit = np.unique(sample_sheet["Description"].values)[0]
       start = bc_kit.find("(") + len("(")
       end = bc_kit.find(")")
       bc_kit = bc_kit[start:end]
    data=dict()
    for index, row in sample_sheet.iterrows():
        assert(row["Sample_ID"] not in data.keys())
        data[row["Sample_ID"]] = dict({"Sample_Name": row["Sample_Name"],
                                       "Sample_Project": row["Sample_Project"],
                                       "barcode_kits": bc_kit,
                                       "index_id": row["I7_Index_ID"],
                                       "Sample_ID": row["Sample_ID"]})
    print("[green] Barcode kit determined as: {}".format(bc_kit))
    return bc_kit, data

def config_to_basecallcmd(config):
    # base cmd.
    cmd = config['guppy_basecaller']['base_calling_cmd']
    # in - and output folder.
    cmd += " -i {}".format(
        config["info_dict"]["poddir"]
    )
    cmd += " -s fastq"
    barcode = False if config['bc_kit'] == 'no_bc' else True
    cmd += " -c {} ".format(config['info_dict']['model'])
    # try:
    #     kit=config["info_dict"]["kit"]
    #     if config["custom_cfg"]:
    #         cmd +=" -c "
    #         if config["protocol"] in ["dna", "cdna"]:
    #             cmd += config["guppy_basecaller"]["dna_model"]+" "
    #             model = config["guppy_basecaller"]["dna_model"]
    #         else:
    #             cmd += config["guppy_basecaller"]["rna_model"]+" "
    #             model = config["guppy_basecaller"]["rna_model"]
    #     else:
    #         cmd +=" -c "
    #         if config["protocol"] in ["dna", "cdna"]:
    #             cmd += config["guppy_basecaller"]["dna_model"]+" "
    #             model = config["guppy_basecaller"]["dna_model"]
    #         else:
    #             cmd += config["guppy_basecaller"]["rna_model"]+" "
    #             model = config["guppy_basecaller"]["rna_model"]
            #flowcell=config["info_dict"]["flowcell"]
            #cmd += " --flowcell {} --kit {} ".format(flowcell,kit)
            #model = " --flowcell {} --kit {} ".format(flowcell,kit)
    # except:
    #     cmd +=" -c "
    #     if config["protocol"] in ["dna", "cdna"]:
    #         cmd += config["guppy_basecaller"]["dna_model"]+ " "
    #         model = config["guppy_basecaller"]["dna_model"]
    #     else:
    #         cmd += config["guppy_basecaller"]["rna_model"]+ " "
    #         model = config["guppy_basecaller"]["rna_model"]
    cmd += config["guppy_basecaller"]["base_calling_options"]
    if barcode is True:
        cmd += " --barcode_kits {} ".format(config["bc_kit"])
        cmd += config["guppy_basecaller"]["base_calling_barcode_options"]
    else:
        cmd += " --trim_strategy dna "
    return (cmd)
