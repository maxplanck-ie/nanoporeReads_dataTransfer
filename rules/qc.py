#!/usr/bin/env python3
import subprocess as sp
import sys
import os
import warnings


rule pycoQc_fastq:
    input:
        rename = "{sample_id}_renamed.done"
    output:
        fastqc = "{sample_id}_qc.done"
    log:
        out ="LOG/fastqc_{sample_id}.log.out",
        err = "LOG/fastqc_{sample_id}.log.err"
    run:
        this_sample = config["data"][wildcards.sample_id]
        if not os.path.exists("FASTQC_Project_"+this_sample["Sample_Project"]):
            os.mkdir("FASTQC_Project_"+this_sample["Sample_Project"])
        fastqc_path = os.path.join(os.path.join("FASTQC_Project_"+this_sample["Sample_Project"],
                                                "Sample_"+this_sample["Sample_ID"]))
        if not os.path.exists(fastqc_path):
            os.mkdir(fastqc_path)
        else:
            warnings.warn("This directory already exists, computing qc is skipped "
                          "and the folder stays untouched.")
            sp.check_call("touch "+output.fastqc, shell=True)
            return

        if this_sample["barcode_kits"] != "no_bc":
            bs = this_sample["index_id"]
            bs = [s for s in " ".join(list(bs)).split() if s.isdigit()]
            bs = "".join(bs)
            print("bs:", bs)
            barcode = "barcode"+bs
            if os.path.isfile("fastq/sequencing_summary_"+barcode+".txt"):
                cmd = config["pycoQc"]["pycoQc"]+" "
                cmd += "fastq/sequencing_summary_"+barcode+".txt "
                cmd += "-o "+fastqc_path+"/pycoqc_"+barcode+".html"
                cmd += " >> "+log.out+" 2> "+log.err
                with open(log.out, "w") as called_cmd:
                    called_cmd.write(cmd)
                called_cmd.close()
                sp.check_call(cmd, shell=True)
                sp.check_call("touch "+output.fastqc, shell=True)
            else:
                with open(log.out, "w") as called_cmd:
                    called_cmd.write("No sequence has been demuxed.")
                called_cmd.close()
                sp.check_call("touch "+output.fastqc, shell=True)

        else:
            assert(len(config["data"])==1)
            cmd = config["pycoQc"]["pycoQc"]+" "
            cmd += "fastq/sequencing_summary.txt "
            cmd += "-o "+fastqc_path+"/pycoqc.html"
            cmd += " >> "+log.out+" 2> "+log.err
            with open(log.out, "w") as called_cmd:
                called_cmd.write(cmd)
            called_cmd.close()

            sp.check_call(cmd, shell=True)
            sp.check_call("touch "+output.fastqc, shell=True)


rule pycoQc_bam:
    input:
        mapped = "{sample_id}.mapped"
    output:
        bam_qc = "{sample_id}.bam.qc"
    log:
        out ="LOG/bamqc_{sample_id}.log.out",
        err = "LOG/bamqc_{sample_id}.log.err"
    run:
        this_sample = config["data"][wildcards.sample_id]
        sample_project = this_sample["Sample_Project"]
        bam = this_sample["Sample_Name"]+".bam"
        group=sample_project.split("_")[-1].lower()
        group_path = os.path.join(config["paths"]["groupDir"],group,"sequencing_data/OxfordNanopore")
        if not os.path.exists(group_path):
            group_path = os.path.join(config["paths"]["external_groupDir"],group,"sequencing_data/OxfordNanopore")
        final_path = os.path.join(group_path, config["input"]["name"])
        # final_path = os.path.join(config["paths"]["groupDir"], group, "sequencing_data/OxfordNanopore",config["input"]["name"])
        path_to_bam = os.path.join(final_path,"Analysis_"+sample_project,"mapping_on_"+config["organism"],"Sample_"+wildcards.sample_id)
        bam = os.path.join(path_to_bam, bam)

        cmd = config["pycoQc"]["pycoQc"]+" "
        if config["bc_kit"]  == 'no_bc':
            cmd += "fastq/sequencing_summary.txt "
        else:
            bs = this_sample["index_id"]
            bs = [s for s in " ".join(list(bs)).split() if s.isdigit()]
            bs = "".join(bs)
            barcode = "barcode"+bs
            if os.path.isfile("fastq/sequencing_summary_"+barcode+".txt"):
                cmd += "fastq/sequencing_summary_"+barcode+".txt "
            else:
                sp.check_call("touch "+output.bam_qc, shell=True)
                return

        cmd += " -a "+bam
        cmd += " -o "+path_to_bam+"/pycoqc_bam.html"
        cmd += " >> "+log.out+" 2> "+log.err
        with open(log.out, "w") as called_cmd:
            called_cmd.write(cmd)
        called_cmd.close()

        sp.check_call(cmd, shell=True)
        sp.check_call("touch "+output.bam_qc, shell=True)

# rule bam_compare: # TODO add NanoComp --bam on all bam files
#     input:
#
#     output:
#
#     run:
