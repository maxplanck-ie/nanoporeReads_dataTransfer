#!/usr/bin/env python3
import sys
import os
import subprocess as sp


rule basecalling:
    input:
        fast5 = config["info_dict"]["fast5"]
    output:
        touch("demux.done")
    log:
        out ="LOG/guppy.log.out",
        err = "LOG/guppy.log.err"

    run:
        cmd = config["guppy_basecaller"]["base_calling_cmd"]
        fast5 = config["info_dict"]["fast5"]
        try:
            os.mkdir("fastq")
        except:
            warnings.warn("fastq folder exists!!")
            return sp.check_call("touch "+log.out+" "+log.err, shell=True)

        fastq = os.path.join("fastq")
        barcode = True
        if config["bc_kit"] == "no_bc":
            barcode = False
        cmd += " -i {} -s {} ".format(fast5, fastq)
        try:
            kit=config["info_dict"]["kit"]
            if config["custom_cfg"]:
                cmd +=" -c "
                if config["protocol"] in ["dna", "cdna"]:
                    cmd += config["guppy_basecaller"]["dna_model"]+" "
                else:
                    cmd += config["guppy_basecaller"]["rna_model"]+" "
            else:
                flowcell=config["info_dict"]["flowcell"]
                cmd += " --flowcell {} --kit {} ".format(flowcell,kit)
        except:
            cmd +=" -c "
            if config["protocol"] in ["dna", "cdna"]:
                cmd += config["guppy_basecaller"]["dna_model"]+" "
            else:
                cmd += config["guppy_basecaller"]["rna_model"]+" "
        cmd += config["guppy_basecaller"]["base_calling_options"]
        if barcode is True:
            cmd += " --barcode_kits {} ".format(config["bc_kit"])
            cmd += config["guppy_basecaller"]["base_calling_barcode_options"]
        else:
            cmd += " -q 0 " # works only for smaples without barcodes
        #Direct RNA sequencing
        if config["protocol"] == "rna":
            cmd += " "+config["guppy_basecaller"]["base_calling_RNA_options"]
        else:
            cmd += " --trim_strategy dna "

        cmd += " >> "+log.out+" 2> "+log.err
        with open(log.out, "w") as called_cmd:
            called_cmd.write(cmd)
        called_cmd.close()
        sp.check_call(cmd, shell=True)
