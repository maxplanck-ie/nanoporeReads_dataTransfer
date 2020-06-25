#!/usr/bin/env python3
import subprocess as sp
import sys
import os

def fastq_qc(config,data):
    for k, v in data.items():
        fastqc = os.path.join(config["info_dict"]["flowcell_path"]+"/FASTQC_Project_"+v["Sample_Project"])
        if not os.path.exists(config["info_dict"]["flowcell_path"]+"/FASTQC_Project_"+v["Sample_Project"]):
           os.mkdir(config["info_dict"]["flowcell_path"]+"/FASTQC_Project_"+v["Sample_Project"])
        os.mkdir(fastqc+"/Sample_"+v["Sample_ID"])
        fastqc = os.path.join(fastqc+"/Sample_"+v["Sample_ID"])
    
        cmd = config["nanocomp"]["qc_cmd"]
        cmd += " --fastq {}".format(config["info_dict"]["flowcell_path"]+"/Project_"+v["Sample_Project"]+"/Sample_"+v["Sample_ID"]+"/"+v["Sample_Name"]+".fastq.gz")
        cmd += " -o {} ".format(fastqc)
        cmd += config["nanocomp"]["qc_options"]+" "+v["Sample_Name"]
        print(cmd)
        sp.check_call(cmd, shell=True)
    return None
