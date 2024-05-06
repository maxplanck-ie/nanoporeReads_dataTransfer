# ont_pipeline.py
# Helper funtions for the ONT pipeline
# (C) 2024 Bioinformatics Core
# Max Plank Institute for Immunobiology and Epigenetics

import sys
import os
import re
import glob
import yaml
import json
import subprocess as sp
import pandas as pd
from rich import print
from npr.communication import send_email
from npr.snakehelper import glob2reports, get_seqdir, guppy2dorado

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
    if os.path.exists(loc1):
        return True
    
    for old_dir in config["paths"]["old_outputDirs"]:
        loc2 = os.path.join(
            old_dir,
            os.path.basename(flowcell),
            'analysis.done'
            )
        if os.path.exists(loc2):
            return True
        
    return False

def filter_flowcell(json, config):
    """
    Decide whether to keep a flowcell based on exclusion rules defined in config
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
    fc_id = fc_dir.split('_')
    if len(fc_id)>1:
        fc_id = fc_id[-2]

    if fc_id in config['ignore']['flowcells']:
        print("ignore fc_id {} because of config ".format(fc_id))
        return True

    # flowcell does not match MPI-IE naming convention: don't filter
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

    # filter flowcells based on config['target_flowcell']
    if  config['target_flowcell']:
        dirs = [d for d in dirs if config['target_flowcell'] in d]

    # Iterate over all flowcell in dir
    for flowcell in dirs:
        if not analysis_done(flowcell, config):
            print(' {} found'.format(flowcell))

        # needed here to communicate flowcell with send_email
        config['info_dict']['base_path'] = flowcell

        # exit if sampleSheet.csv does not exists
        ss = os.path.join(flowcell, 'SampleSheet.csv')

        ss2 = os.path.join(
                    config['paths']['outputDir'],
                    os.path.basename(flowcell),
                    'reports',
                    'SampleSheet.csv')

        if not os.path.isfile(ss):
            # here we need some method to retrieve the SampleSheet.csv from Parkour
            # get_samplesheet_from_parkour(flowcell, config)
            if not os.path.isfile(ss2):
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
     - add flags for basecalling and modbed
    """

    # adding default flags for basecalling, align and modbed
    info_dict["do_basecall"] = config["default_process"]["do_basecall"]
    info_dict["do_align"] = config["default_process"]["do_align"]
    info_dict["do_modbed"] = config["default_process"]["do_modbed"]
    
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
    html_file = glob.glob(
        os.path.join(
            flowcell_path,
            'reports',
            '*html'
        )
    )
    if html_file:
        match_html_key_val = re.compile(r'"title": "(.+?)", "value": "(.+?)"')

        print("[green]Reading info from html[/green]")
        html_cont = ''
        with open(html_file[0]) as f:
            html_cont = f.read()
        
        info_dict["software"] = {}
        matches = re.finditer(match_html_key_val, html_cont)
        for match in matches:
            name = match.group(1)
            value = match.group(2)
            if name == "Flow cell type":
                info_dict["flowcell"] = value
            elif name == "Kit type":
                info_dict["kit"] = value
            elif name == "MinKNOW":
                info_dict["software"]["MinKNOW"] = value
            elif name == "Bream":
                info_dict["software"]["Bream"] = value
            elif name == "MinKNOW Core":
                info_dict["software"]["MinKNOW Core"] = value
            elif name == "Configuration":
                info_dict["software"]["Configuration"] = value
            elif name == "Dorado":
                info_dict["software"]["Dorado"] = value

    if json_file:
        print("[green]Reading info from json[/green]")
        with open(json_file[0]) as f: # assume only 1 file
            jsondata = json.load(f)

        if not "flowcell" in info_dict:
            info_dict["flowcell"] = jsondata['protocol_run_info']['meta_info']['tags']['flow cell']['string_value']
        if not "kit" in info_dict:
            info_dict["kit"] = jsondata['protocol_run_info']['meta_info']['tags']['kit']['string_value']
        # HTML file is not reporting barcoding, we need it from the json
        info_dict['barcoding'] = bool(jsondata['protocol_run_info']['meta_info']['tags']['barcoding']['bool_value'])
        
        # check run parameters to capture model, check base calling and alignment
        # model is better defined in json, in html is badly reported
        model = None
        
        for par in jsondata['protocol_run_info']['args']:
            if par.startswith('--model_filename='):
                px,model = par.split("=")
                info_dict['model_def'] = model
                info_dict['model'] = model
                print (f"  [green]Found model as {model}[/green]")
            elif par == '--base_calling=on':
                info_dict['do_basecall'] = 'no_basecall'
                print (f"  [green]Found basecalling is already done[/green]")
            elif par == '--alignment':
                info_dict['do_align'] = 'no_align'
                print (f"  [green]Found alignment is already done[/green]")
            
        if not model:
            print("Model was not found in command parameters, capturing the default value")      
            info_dict['model_def'] = jsondata['protocol_run_info']['meta_info']['tags']['default basecall model']['string_value']
            info_dict['model'] = info_dict['model_def']

        if 'modbases' in info_dict['model_def']:
            info_dict['do_modbed'] = 'do_modbed'
            print (f"  [green]Found modified bases was used, BED will be genered[/green]")
        

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

        # getting software and versions unless already defined by the html
        if not "software" in info_dict:
            if 'software_versions' in jsondata['protocol_run_info']:
                info_dict['software'] = jsondata['protocol_run_info']['software_versions']
                if 'guppy_connected_version' in info_dict['software']:
                    info_dict['software']['Dorado'] = info_dict['software']['guppy_connected_version']

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

    if info_dict["do_basecall"] == "do_basecall":
        if (config['basecaller']=="guppy"):
            if config['guppy_basecaller']['guppy_model'] is not None:
                # record _full_ absolute path to model
                info_dict['model'] = config['guppy_basecaller']['guppy_model']
            else:
                # infer model path name from flowcell and kit information using dictionary
                modeldic = yaml.safe_load(
                    open(config['guppy_basecaller']['model_dictionary'])
                )
                info_dict['model'] = modeldic[info_dict['flowcell']][info_dict['kit']]
                # modify model name if modification calling is desired
                if config['guppy_basecaller']['guppy_mod'] is not None:
                    patt = r'(\w*)_(\w{3,4}\.cfg)'
                    repl = "modbases_" + config['guppy_basecaller']['guppy_mod']
                    info_dict['model'] = re.sub(patt, r'\1_{}_\2'.format(repl),info_dict['model'])

        elif (config['basecaller']=="dorado"):
            if config['dorado_basecaller']['dorado_model'] is not None:
                # record _full_ abnsolute path to model
                info_dict['model'] = config['dorado_basecaller']['dorado_model']
            else:
                # default name of model derived from json (see above)
                model_name = info_dict['model_def']
                info_dict['model'] = os.path.join(
                    config['dorado_basecaller']['model_directory'],
                    model_name
                )

    print('model = {}'.format(info_dict['model']))

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

    # convert: No_index? --> no_bc
    if 'I7_Index_ID' not in sample_sheet.columns:
        print("[red] Column I7_Index_ID not defined in sample_sheet [/red]")
        exit(1)
    # legacy fix: replace missing barcode (NaN or 'No_index.*') by 'no_bc'
    sample_sheet['I7_Index_ID'] = sample_sheet['I7_Index_ID'].fillna('no_bc')  
    sample_sheet['I7_Index_ID'] = sample_sheet['I7_Index_ID'].replace('No_index.*','no_bc', regex=True)

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
    Warning1: assumes a _single_ project per flowcell
    Warning2: assumes regular group directory - not applicable to 'external' projects
    """
    group    = config['data']['projects'][0].split("_")[-1].lower()
    groupdir = os.path.join(config["paths"]["groupDir"],group)
    # Warning: get_seqdir _creates_ directories as side effect
    groupONT = get_seqdir(groupdir, "sequencing_data")
    fc_base  = os.path.basename(config['info_dict']['flowcell_path'])
    periphery = os.path.join(groupONT,fc_base)
    return(periphery)


def get_dest_path(config, dir):
    '''
    for a given directory 'dir' (Project_\d+)_(\w+)_(\w+))
    get target destination based on 'pi_name' in 'dir'
    in contrast to get_periphery, this function can handle "unknown" or 
    undeterminable PIs and will invoke an external_groupDir in such cases
    '''

    # default pi_name if it cannot be determined from basename(dir)
    pi_name = "unknown"

    # try to get pi_name from basename(dir)
    match = re.match(r'^(Project_\d+)_(\w+)_(\w+)', os.path.basename(dir))
    if match:
        p_id, username, pi_name = match.groups()
        pi_name = pi_name.lower()

    # assume PI has directories all set up
    groupdir = os.path.join(config["paths"]["groupDir"],pi_name)
    if not os.path.exists(groupdir):
        # group pi_name has no volumne so transfer somewhere else (external runs)
        dest_path = os.path.join(
            config["paths"]["external_groupDir"],
            pi_name,
            "sequencing_data/OxfordNanopore"
        )
        print(f"No group directory {groupdir} for PI {pi_name}")
        print(f"Use {dest_path}")
    else:
        # get the most recent destination path: e.g. "{groupdir}/sequencing_data6/OxfordNanopore"
        dest_path = get_seqdir(groupdir, "sequencing_data")

    dest_path = os.path.join(
        dest_path,
        os.path.basename(config['info_dict']['flowcell_path'])
    )
    return dest_path
