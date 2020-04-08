#!/usr/bin/env python3
import sys
import os
import subprocess as sp

def base_calling(config):
    """
    run guppy basecaller
    """
    cmd = config["guppy_basecaller"]["base_calling_cmd"]
    fast5 = config["info_dict"]["fast5"]
    os.mkdir(config["info_dict"]["flowcell_path"]+"/Project_"+config["data"]["Sample_Project"])
    fastq = os.path.join(config["info_dict"]["flowcell_path"],"Project_"+config["data"]["Sample_Project"])
    os.mkdir(fastq+"/Sample_"+config["data"]["Sample_Name"])
    fastq = os.path.join(fastq,"Sample_"+config["data"]["Sample_Name"])
    config["info_dict"]["fastq"] = fastq
    flowcell=config["info_dict"]["flowcell"]
    barcode = False
    if config["data"]["barcode_kits"] is not "no_bc": #At the moment there is no multi samples per run so i dunno how would fast5 look and how i should update the code for that
        barcode = True
    kit=config["info_dict"]["kit"]
    cmd += " -i {} -s {} --flowcell {} --kit {} ".format(fast5,fastq,flowcell,kit)
    cmd += config["guppy_basecaller"]["base_calling_options"]
    if barcode is True:
        cmd += " --barcode_kits {} ".format(kit)
        cmd += config["guppy_basecaller"]["base_calling_barcode_options"]
    else:
        cmd += " -q 0 " # TODO test on sampels with barcodes
    # RNA
    if "RNA" in config["info_dict"]["kit"]:
        cmd += " "+config["guppy_basecaller"]["base_calling_RNA_options"]
    sp.check_call(cmd, shell=True)
