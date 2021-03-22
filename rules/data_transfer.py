#!/usr/bin/env python3

import subprocess as sp
import sys
import os
import warnings
import shutil


def genome_index(config, path):
    """
    Generating genome indices for minimap2
    """
    reference = config["organism"]
    reference_fa = config["genome"][reference]
    cmd = "ln -s " + reference_fa + " " + path + "/" + reference + "_genome.fa;"
    cmd += config["mapping"]["mapping_cmd"] +" "+ config["mapping"]["index_options"] + " "
    cmd += path + "/" + reference + "_genome.mmi "
    cmd += path + "/" + reference + "_genome.fa "
    print(cmd)
    sp.check_call(cmd, shell=True)


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
            genome_index(config, analysis+"/mapping_on_"+config["organism"])
        analysis_dir = analysis+"/mapping_on_"+config["organism"]
        os.mkdir(analysis_dir+"/Sample_"+wildcards.sample_id)
        sp.check_call("touch "+output.transferred, shell = True)
