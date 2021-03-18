#!/usr/bin/env python3

import subprocess as sp
import sys
import os
import warnings
import shutil


# Trnasfer Project_ and FASTQC_Project_ to the user's directory
rule data_transfer:
    input:
        fastqc = "{sample_id}_qc.done"
    output:
        transferred = temp("{sample_id}.transferred")
    run:
        this_sample = config["data"][wildcards.sample_id]
        group=this_sample["Sample_Project"].split("_")[-1].lower()
        final_path = os.path.join(config["paths"]["groupDir"],group,"sequencing_data/OxfordNanopore/"+config["input"]["name"])
        if not os.path.exists(final_path):
            os.mkdir(final_path)
        if not os.path.exists(final_path+"/Project_"+this_sample["Sample_Project"]):
            fastq = config["info_dict"]["flowcell_path"]+"/Project_"+this_sample["Sample_Project"]
            shutil.copytree(fastq,final_path+"/Project_"+this_sample["Sample_Project"])
        if not os.path.exists(final_path+"/FASTQC_Project_"+this_sample["Sample_Project"]):
            fastqc = config["info_dict"]["flowcell_path"]+"/FASTQC_Project_"+this_sample["Sample_Project"]
            shutil.copytree(fastqc,final_path+"/FASTQC_Project_"+this_sample["Sample_Project"])
        analysis = os.path.join(final_path,"Analysis_"+this_sample["Sample_Project"])
        if not os.path.exists(analysis):
            os.mkdir(analysis)
            os.mkdir(analysis+"/mapping_on_"+config["organism"])
        analysis = analysis+"/mapping_on_"+config["organism"]
        os.mkdir(analysis+"/Sample_"+wildcards.sample_id)
        sp.check_call("touch "+output.transferred, shell = True)
