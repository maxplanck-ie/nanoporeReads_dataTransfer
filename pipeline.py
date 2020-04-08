#!/usr/bin/env python3
import sys
import os
import argparse
import configparser
import shutil
import pandas as pd
import subprocess as sp
from qc import *
from mapping import *
from baseCalling import base_calling

def get_parser():

    parser = argparse.ArgumentParser(description='A Pipeline to process fast5.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # required argumnets:
    required.add_argument("-i",
                        "--input",
                        type=str,
                        dest="input",
                        help='input path')
    required.add_argument("-r",
                        "--ref",
                        type=str,
                        dest="reference",
                        help='reference genome')
    return parser

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
        sys.exit("a flowcell with the same ID already exists!!")
    flowcell_path = os.path.join(config["paths"]["outputDir"]+input)
    info_dict["flowcell_path"] = flowcell_path
    if not os.path.exists(flowcell_path+"/fast5"):
         sys.exit("fast5 path doesnt exist.")
    info_dict["fast5"] = os.path.join(flowcell_path,"fast5")

    summary_file = [filename for filename in os.listdir(flowcell_path) if filename.startswith("final_summary")]
    if summary_file is []:
         sys.exit("final summary file doesnt exist.")
    assert len(summary_file) == 1
    summary_file = os.path.join(flowcell_path,summary_file[0])
    with open(summary_file,"r") as f:
        for line in f.readlines():
            if line.startswith("protocol="):
                info_dict["flowcell"] = line.split(":")[1]
                if info_dict["flowcell"] not in config["flowcell"]["compatible_flowcells"]:
                    sys.exit("flowcell id is not valid!")
                info_dict["kit"] = line.split(":")[2]
                if info_dict["kit"] not in config["flowcell"]["compatible_kits"]:
                    sys.exit("kit id is not valid!")
                return info_dict

    return None

def read_samplesheet(config, ref_genome):
    """
        read samplesheet
    """
    sample_sheet = pd.read_csv(config["info_dict"]["flowcell_path"]+"/Samplesheet.csv",
                               sep = ",", skiprows=[0])
    sample_sheet = sample_sheet.fillna("no_bc")
    config["data"]=dict()
    for index, row in sample_sheet.iterrows():
        if index == 0:
            config["data"]["Sample_Name"] = row["Sample_Name"]
            config["data"]["Sample_Project"] = row["Sample_Project"]
            config["data"]["barcode_kits"] = row["barcode_kits"]
        else:
            config["data"]["Sample_Name"].append(row["Sample_Name"])
            config["data"]["Sample_Project"].append(row["Sample_Project"])
            config["data"]["barcode_kits"].append(row["barcode_kits"])
    config["data"]["ref"] = ref_genome

def rename_fastq(config):
    fastq_path = config["info_dict"]["fastq"]
    sample_name = config["data"]["Sample_Name"]
    if config["data"]["barcode_kits"] is not "no_bc":
        barcode = config["data"]["barcode_kits"]
        fastq_path = os.path.join(fastq_path,barcode)
        config["info_dict"]["fastq"] = fastq_path

    cmd = "cat {}/*.fastq.gz > {}/{}.fastq.gz ;".format(fastq_path,fastq_path,sample_name)
    cmd += "rm {}/fastq_runid*.fastq.gz".format(fastq_path)
    sp.check_call(cmd, shell=True)


def transfer_data(config):
    """
    Trnasfer Project_ and FASTQC_Project_ to the user's directory
    """
    group=config["data"]["Sample_Project"].split("_")[2]
    final_path = "/"+group+"/sequencing_data/"+config["input"]["name"]
    if not os.path.exists(config["paths"]["groupDir"]+final_path):
        os.mkdir(config["paths"]["groupDir"]+final_path)
        final_path = os.path.join(config["paths"]["groupDir"],final_path)
    else:
        sys.exit("a flowcell with the same ID already exists!!")
    fastq = config["info_dict"]["flowcell_path"]+"/Project_"+config["data"]["Sample_Project"]
    fastqc = config["info_dict"]["flowcell_path"]+"/FASTQC_Project_"+config["data"]["Sample_Project"]
    shutil.copytree(fastq,final_path+"/Project_"+config["data"]["Sample_Project"])
    shutil.copytree(fastqc,final_path+"/FASTQC_Project_"+config["data"]["Sample_Project"])
    analysis = final_path+"/Analysis_"+config["data"]["Sample_Project"]
    os.mkdir(analysis)
    os.mkdir(analysis+"/mapping_on_"+config["data"]["ref"])
    config["data"]["analysis"] = analysis
    config["data"]["mapping"] = analysis+"/mapping_on_"+config["data"]["ref"]



def main():
    args = get_parser().parse_args()

    config = configparser.ConfigParser()
    config.read_file(open("./config.ini","r"))
    config["input"]=dict([("name",os.path.basename(os.path.realpath(args.input)))])
    info_dict = read_flowcell_info(config)
    config["info_dict"]=info_dict
    read_samplesheet(config, args.reference)
    base_calling(config)
    rename_fastq(config)
    fastq_qc(config)
    transfer_data(config)
    if "RNA" in config["info_dict"]["kit"]:
        mapping_rna(config)
    else:
        mapping_dna(config)


if __name__== "__main__":
    main()
