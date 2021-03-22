#!/usr/bin/env python3
import sys
import os
import subprocess as sp


rule basecalling:
    input:
        fast5 = config["info_dict"]["fast5"]
    output:
        fastq = directory("fastq"),
        demux = temp("demux.done")
    run:
        cmd = config["guppy_basecaller"]["base_calling_cmd"]
        fast5 = config["info_dict"]["fast5"]
        os.mkdir(output.fastq)
        fastq = os.path.join(output.fastq)
        flowcell=config["info_dict"]["flowcell"]
        barcode = True
        if config["bc_kit"] == "no_bc":
            barcode = False
        kit=config["info_dict"]["kit"]
        cmd += " -i {} -s {} --flowcell {} --kit {} ".format(fast5,
                                                             fastq,
                                                             flowcell,kit)
        cmd += config["guppy_basecaller"]["base_calling_options"]
        if barcode is True:
            cmd += " --barcode_kits {} ".format(config["bc_kit"])
            cmd += config["guppy_basecaller"]["base_calling_barcode_options"]
        else:
            cmd += " -q 0 " # works only for smaples without barcodes
        #Direct RNA sequencing
        if "RNA" in config["info_dict"]["kit"]:
            cmd += " "+config["guppy_basecaller"]["base_calling_RNA_options"]

        print(cmd)
        sp.check_call(cmd, shell=True)
        sp.check_call("touch "+output.demux, shell=True)
