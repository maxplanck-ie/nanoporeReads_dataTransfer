#!/usr/bin/env python3
import subprocess as sp
import sys
import os

def fastq_qc(config):
    sample_name = config["data"]["Sample_Name"]
    os.mkdir(config["info_dict"]["flowcell_path"]+"/FASTQC_Project_"+config["data"]["Sample_Project"])
    fastqc = os.path.join(config["info_dict"]["flowcell_path"],"FASTQC_Project_"+config["data"]["Sample_Project"])
    os.mkdir(fastqc+"/Sample_"+config["data"]["Sample_Name"])
    fastqc = os.path.join(fastqc,"Sample_"+config["data"]["Sample_Name"])
    config["info_dict"]["fastqc"] = fastqc
    cmd = config["nanocomp"]["qc_cmd"]
    cmd += " --fastq {}".format(config["info_dict"]["fastq"]+"/"+sample_name+".fastq.gz")
    cmd += " -o {} ".format(fastqc)
    cmd += config["nanocomp"]["qc_options"]+" "+sample_name
    sp.check_call(cmd, shell=True)
    return None
