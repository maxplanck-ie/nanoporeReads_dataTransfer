import sys
import os
import re
import glob
import pandas as pd
import numpy as np
from rich import print
from importlib.metadata import version
import yaml
import json
import subprocess as sp
from npr.communication import send_email
from npr.snakehelper import glob2reports, get_seqdir

def analysis_done(flowcell, config):
    """
    Determine whether flowcell has already been analysed
    Notice that such an analysis could have been done at different locations
    """
    loc1 = os.path.join(
        config["paths"]["outputDir"],
        os.path.basename(flowcell),
        'analysis.done'
        )
    loc2 = os.path.join(
        config["paths"]["old_outputDir"],
        os.path.basename(flowcell),
        'analysis.done'
        )

    if os.path.exists(loc1) or os.path.exists(loc2):
        return True

    return False

def filter_flowcell(json, config):
    """
    decide whether to keep a flowcell based on exclusion rules defined in config
    the exclusion rules are specific to our MPI-IE setup

    example json file: "report_PAQ97481_20230818_1514_e1253480.json"
    example fc direct: "20230818_1512_P2S-00500-A_PAQ97481_e1253480"
    """

    # get parent directory of json report = flowcell directory
    fc_dir = os.path.dirname(json)
    fc_dir = os.path.basename(fc_dir.rstrip('/'))
    #print("json: {} , fc_dir: {} ".format(json, fc_dir))

    if fc_dir in config['ignore']['dirs']:
        print("ignore fc_dir {} because of config ".format(fc_dir))
        return True
    
    # MPI-IE specific pattern for flowcells
    pattern = r"\d{8}_\d{4}_\w+-\d{5}-[A-Z0-9]+_([A-Z0-9]+)_\d+"
    match = re.match(pattern, fc_dir)
    if match:
        fc_id = match.group(1)
        if fc_id in config['ignore']['flowcells']:
            print("ignore fc_id {} because of config ".format(fc_id))
            return True
        
    # flowcell does not match MPI-IE naming convention: don't filter
    return False

    # be extra cautious: get fc_id and fc_dir from (large) json file
    # this is slow, but it will only be done for a few non-standard flowcells
    fc_yaml = yaml.safe_load(open(json))
    fc_id  = fc_yaml['protocol_run_info']['flow_cell']['flow_cell_id']
    fc_dir = fc_yaml['protocol_run_info']['output_path']
    fc_dir = os.path.basename(fc_dir.rstrip('/'))  
    if fc_dir in config['ignore']['dirs']:
        print("ignore fc_dir {} because of config ".format(fc_dir))
        return True
    if fc_id in config['ignore']['flowcells']:
        print("ignore fc_id {} because of config ".format(fc_id))
        return True

    return False


def find_new_flowcell(config):
    """
    look for new flowcells inf offloadDir   
    baseDir flowcells are flagged with 'transfer.final'
    """
    #################################### offload ####################################
    offload_path = os.path.join(
        config['paths']['offloadDir']
    )
    dirs = []
    if offload_path:
        # Assume that report*.json mark 'ready' flow cells
        # depending on the offload_path
        pattern = os.path.join(offload_path, '**', 'report*.json')
        jsons = glob.glob(pattern, recursive=True)

        for j in jsons:
            if not filter_flowcell(j,config):
                # collect full path to json
                dirs.append(os.path.dirname(j))

    # Iterate over all flowcell in dir
    for flowcell in dirs:
        print('Working with {}'.format(flowcell))

        # test and continue if 'analysis.done' exists for this flowcell
        if analysis_done(flowcell, config):
            continue

        # exit if sampleSheet.csv does not exists
        ss = os.path.join(flowcell, 'SampleSheet.csv')
        if not os.path.isfile(ss):
            msg = "No SampleSheet.csv file.\n"
            send_email('Error for flowcell:', msg, config)
            sys.exit("no sampleSheet.")
        
        # return flowcell to ont()
        msg = "SampleSheet.csv file found.\n"
        config["input"] = {
            'name': os.path.basename(flowcell)
        }
        return (os.path.basename(flowcell), msg, flowcell)

    return (None, None, None)

def read_flowcell_info(config, info_dict, base_path):
    """
     - set paths
     - copy json, summaries & samplesheet
     - parse json for flowcell, kit & barcoding
     - infer model from flowcell + kit
    """
    flowcell = config["input"]["name"]
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
        print("[green]Reading info from json[/green]")
        with open(json_file[0]) as f: # assume only 1 file
            jsondata = json.load(f)
        info_dict["flowcell"] = jsondata['protocol_run_info']['meta_info']['tags']['flow cell']['string_value']
        info_dict["kit"] = jsondata['protocol_run_info']['meta_info']['tags']['kit']['string_value']
        info_dict['barcoding'] = bool(jsondata['protocol_run_info']['meta_info']['tags']['barcoding']['bool_value'])
        info_dict['model_def'] = jsondata['protocol_run_info']['meta_info']['tags']['default basecall model']['string_value']

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
            if info_dict['barcoding']:
                print('[red] Not barcoding kit found in json. Default to flowcell kit.[/red]')
                info_dict['barcode_kit'] = info_dict['kit']
            else:
                info_dict['barcode_kit'] = 'no_bc'
    else:
        # try to get txt file.
        print('base path == {}'.format(base_path))
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
            seq_sum = glob.glob(
                os.path.join(
                    base_path,
                    'sequencing_summary*txt'
                )
            )
            if not seq_sum:
                sys.exit("No sequencing summary.txt file found. exiting.")
            else:
                # We try and fetch a barcode kit from sequencing summary.
                # only open up the first 2 lines, as these can grow long.
                # if barcoding is there, there'll be a barcoding
                lines = 0
                head = []
                with open(seq_sum[0], 'r') as f:
                    for line in f:
                        if lines < 2:
                            head.append(line.strip().split())
                            lines += 1
                        else:
                            break
                    if 'barcode_kit' in head[0]:
                        info_dict['barcoding'] = True
                        bkit = head[1][head[0].index('barcode_kit')]
                        if info_dict['kit'] == 'SQK-PCB111-24' and info_dict['barcoding']:
                            info_dict['barcode_kit'] = 'SQK-PCB111-24'
                        else:
                            info_dict['barcode_kit'] = bkit
                        print('Barcoding detected: kit = {}'.format(bkit))
                    else:
                        print('no evidence for barcoding in seq summary.')
                        info_dict['barcode_kit'] = 'no_bc'
                        info_dict['barcoding'] = False
                    print(head[0].index('barcode_kit'))
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
    if len(data['samples']) > 1 and config["info_dict"]['barcoding'] == False:
        print("[red] Danger, barcoding inferred as false, but more then 1 sample in samplesheet. [/red]")
        print("[red] Overwrite barcoding boolean to true. [/red]")
        config["info_dict"]['barcoding'] = True
        print("[red] Defaulting to flowcell kit as barcode kit ! [/red]")
        config["info_dict"]['barcode_kit'] = config["info_dict"]['kit']
    print("[green] Barcode kit determined as: {}".format(
        config["info_dict"]['barcode_kit']
    ))
    return config["info_dict"]['barcode_kit'], data

def get_periphery(config):
    """
    Return the full pathname to the flowcell in the periphery
    as it should be there after transfer
    """
    group    = config['data']['projects'][0].split("_")[-1].lower()
    groupdir = os.path.join(config["paths"]["groupDir"],group)
    # Warning: get_seqdir _creates_ directories as side effect
    groupONT = get_seqdir(groupdir, "sequencing_data")
    fc_base  = os.path.basename(config['info_dict']['flowcell_path'])
    periphery = os.path.join(groupONT,fc_base)
    return(periphery)
