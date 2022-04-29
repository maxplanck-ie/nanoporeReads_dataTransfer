#!/usr/bin/env python3

import subprocess as sp
import sys
import os
import warnings
import glob
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
def get_seqdir(groupdir, seqdir):
    max_num = 0
    for dir in glob.glob(os.path.join(groupdir, seqdir+"*")):
        if dir.split(seqdir)[-1] != "":
            num = dir.split(seqdir)[-1]
            max_num = max(max_num, int(num))
    if max_num != 0:
        seqdir = seqdir + str(max_num)
    if not os.path.exists(os.path.join(groupdir, seqdir, "OxfordNanopore")):
        os.makedirs(os.path.join(groupdir, seqdir, "OxfordNanopore"))
    return os.path.join(groupdir, seqdir, "OxfordNanopore")

# Trnasfer Project_ and FASTQC_Project_ to the user's directory and generates the Analysis_ folder directly under the users path
rule data_transfer:
    input:
        fastqc = "{sample_id}_qc.done"
    output:
        transferred = touch("{sample_id}.transferred")
    run:
        this_sample = config["data"][wildcards.sample_id]
        group=this_sample["Sample_Project"].split("_")[-1].lower()
        groupdir = os.path.join(config["paths"]["groupDir"],group)
        if not os.path.exists(groupdir): # If external (no /data/pi/ path exists)
            group_path = os.path.join(config["paths"]["external_groupDir"],group,
                                      "sequencing_data/OxfordNanopore")
            if not os.path.exists(group_path):
                os.makedirs(group_path)
        else:
            group_path = get_seqdir(groupdir, "sequencing_data")
        final_path = os.path.join(group_path, config["input"]["name"])
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
        if not os.path.exists(os.path.join(analysis, "mapping_on_"+config["organism"])):
            os.mkdir(os.path.join(analysis, "mapping_on_"+config["organism"]))
        genome_index(config, os.path.join(analysis, "mapping_on_"+config["organism"]))
        analysis_dir = os.path.join(analysis,"mapping_on_"+config["organism"])
        try:
            os.mkdir(os.path.join(analysis_dir,"Sample_"+wildcards.sample_id))
        except:
            warnings.warn("{} already exists!!".format(os.path.join(analysis_dir,"Sample_"+wildcards.sample_id)))
