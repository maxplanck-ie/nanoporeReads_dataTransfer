#!/usr/bin/env python3
import sys
import os
import warnings
import argparse
import configparser
import shutil
import pandas as pd
import subprocess as sp
import yaml
import glob
from threading import Event
import signal
import requests


gotHUP = Event()

def breakSleep(signo, _frame):
    gotHUP.set()


def sleep(config):
    print("go to sleep")
    gotHUP.wait(timeout=float(config['options']['sleep_time'])*60*60) # in second
    gotHUP.clear()

signal.signal(signal.SIGHUP, breakSleep)


def get_parser(): # TODO can be removed!

    parser = argparse.ArgumentParser(description='A Pipeline to process fast5.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    return parser


def find_new_flowcell(config):
    """
    look for new flowcells with `transfer.final` and `SampleSheet.csv`
    """
    base_path = os.path.join(config["paths"]["baseDir"])
    dirs = glob.glob(os.path.join(base_path, "*/transfer.final"))
    for dir in dirs:
        flowcell = os.path.dirname(dir)
        # first make sure the dir is not in the baseDir.
        if not os.path.basename(flowcell) in os.listdir(config["paths"]["outputDir"]):
            print("New flowcell found: {}".format(flowcell))
            if not os.path.isfile(os.path.join(flowcell, 'SampleSheet.csv')):
                print('No sampleSheet for {}. exiting.'.format(flowcell))
                sys.exit()
            config["input"]=dict([("name",os.path.basename(flowcell))])
            return(os.path.basename(flowcell))
    sleep(config)


#        if os.path.isfile(os.path.join(flowcell, 'SampleSheet.csv')):
#            if os.path.basename(flowcell) in os.listdir(config["paths"]["outputDir"]):
#                print("flowcell already exists")
#                continue
#            else:
#                print("I found a new flowcell!")
#                config["input"]=dict([("name",os.path.basename(flowcell))])
#                return os.path.basename(flowcell)
#        else:
#            print("there is no samplesheet")
#            sleep(config)
#    sleep(config)

def query_parkour(config, flowcell):
    """

    """
    fc = flowcell.split("_")[3]
    d = {'flowcell_id': fc}
    res = requests.get(config["parkour"]["url"],\
                       auth=(config["parkour"]["user"], config["parkour"]["password"]), params=d)
    if res.status_code == 200: # Flowcell exists!
        info_dict = res.json()
        print(info_dict)
        first_key = list(info_dict.keys())[0]
        first_entry = list(info_dict[first_key].keys())[0]
        organism = info_dict[first_key][first_entry][-3]
        protocol = info_dict[first_key][first_entry][-2]
        if 'cDNA' in protocol:
            protocol = 'cdna'
        elif 'DNA' in protocol:
            protocol = 'dna'
        elif "RNA" in protocol:
            protocol = 'RNA'
        else:
            protocol = 'dna' # TODO just to let the pipeline run for now!
            # exit("protocol not found")
        config["organism"] = str(organism)
        config["protocol"] = protocol
        print(protocol)
    else:
        print("flowcell does not exist")



def read_flowcell_info(config):
    """
    Check the flowcell path to find the info needed for base calling
    """
    input = config["input"]["name"]
    info_dict = dict()
    base_path = os.path.join(config["paths"]["baseDir"]+input)
    if not os.path.exists(config["paths"]["outputDir"]+input):
        shutil.copytree(base_path,config["paths"]["outputDir"]+input)
    else:
        warnings.warn("a flowcell with the same ID already exists!! Does not overwrite") # TODO should be changed to sys.exit
    flowcell_path = os.path.join(config["paths"]["outputDir"]+input)
    info_dict["flowcell_path"] = flowcell_path
    if glob.glob(flowcell_path+"/fast5*") == []: # Use both fast5_pass and fast5_fail
         sys.exit("fast5 path doesnt exist.")
    info_dict["fast5"] = os.path.join(flowcell_path)

    summary_file = [filename for filename in os.listdir(flowcell_path) if filename.startswith("final_summary")]
    if summary_file == []:
         warnings.warn("final summary file doesnt exist. It uses the -c in guppy.")
         config["no_kit_info"] = True
         config["custom_cfg"] = True
    else:
        config["no_kit_info"] = False
        config["custom_cfg"] = False
        assert len(summary_file) == 1
        summary_file = os.path.join(flowcell_path,summary_file[0])
        with open(summary_file,"r") as f:
            for line in f.readlines():
                if line.startswith("protocol="):
                    try:
                        info_dict["flowcell"] = line.split(":")[1]
                        if info_dict["flowcell"] not in config["flowcell"]["compatible_flowcells"]:
                            sys.exit("flowcell id is not valid!")
                        info_dict["kit"] = line.split(":")[2]
                        if info_dict["kit"].endswith("\n"):
                            info_dict["kit"] = info_dict["kit"].split("\n")[0]
                        if str(info_dict["kit"]) not in config["flowcell"]["compatible_kits"]:
                            sys.exit("kit id is not valid!")
                    except:
                        warnings.warn("final summary format is incorrect! It uses the -c in guppy.")
                        config["no_kit_info"] = True
                        config["custom_cfg"] = True

    return info_dict


def read_samplesheet(config):
    """
        read samplesheet
    """
    sample_sheet = pd.read_csv(config["info_dict"]["flowcell_path"]+"/SampleSheet.csv",
                               sep = ",", skiprows=[0])
    # sample_sheet = sample_sheet.fillna("no_bc")
    sample_sheet['I7_Index_ID'] = sample_sheet['I7_Index_ID'].str.replace('No_index1','no_bc', regex = True) # TODO!! need to be applied on bc kit too!
    # assert(len(sample_sheet["barcode_kits"].unique())==1)
    # bc_kit = sample_sheet["barcode_kits"].unique()[0]

    if any(sample_sheet['I7_Index_ID'].str.contains('no_bc')):
        bc_kit = "no_bc"
    else:
        bc_kit = "SQK-PCB109" # TODO just for testing
    print(sample_sheet)
    data=dict()
    for index, row in sample_sheet.iterrows():
        assert(row["Sample_ID"] not in data.keys())
        data[row["Sample_ID"]] = dict({"Sample_Name": row["Sample_Name"],
                                       "Sample_Project": row["Sample_Project"],
                                       # "barcode_kits": row["barcode_kits"], TODO
                                       "barcode_kits": bc_kit, # TODO just for testing
                                       "index_id": row["I7_Index_ID"],
                                       "Sample_ID": row["Sample_ID"]})
    print(bc_kit)
    return bc_kit, data



def report_contamination(config, data, protocol):
    if protocol == 'rna':
        mapping_rna_contamination(config, data)



def main():
    # parse arguments
    args = get_parser().parse_args()

    # read config
    config = yaml.safe_load(open(os.path.join(os.path.dirname(__file__), 'config.yaml')))
    root = config["paths"]["baseDir"]

    flowcell = find_new_flowcell(config)
    query_parkour(config, flowcell)

    # read the flowcell info & copy it over from dont_touch_this to rapidus
    info_dict = read_flowcell_info(config)
    config["info_dict"]=info_dict

    # read samplesheet
    bc_kit,data = read_samplesheet(config)
    config["data"] = data
    config["bc_kit"] = bc_kit
    print(config["data"])

    # write the updated config file under the output path
    configFile = os.path.join(config["paths"]["outputDir"], config["input"]["name"], "pipeline_config.yaml")
    with open(configFile, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    #run snakemake
    output_directory = os.path.join(config["paths"]["outputDir"], config["input"]["name"])
    snakefile_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), "ont_pipeline.Snakefile")
    # snakemake log file
    fnames = glob.glob(os.path.join(output_directory, 'ont_run-[0-9]*.log'))
    if len(fnames) == 0:
        n = 1  # no matching files, then this is the first run
    else:
        fnames.sort(key=os.path.getctime)
        n = int(fnames[-1].split("-")[-1].split(".")[0]) + 1  # get new run number
    # append the new run number to the file name
    logfile_name = "ont_run-{}.log".format(n)
    snakemake_cmd = " snakemake  -s "+snakefile_directory+" --jobs 5 -p --verbose \
                     --configfile "+configFile+" \
                     --directory " + output_directory  \
                     + " --debug-dag " \
                     +" 2> "+os.path.join(output_directory, logfile_name)

                     # + " --dag | dot -Tpdf > dag.pdf"

    sp.check_output(snakemake_cmd, shell = True)
    sleep(config)

if __name__== "__main__":
    main()
