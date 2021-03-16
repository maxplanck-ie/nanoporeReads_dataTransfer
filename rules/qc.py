#!/usr/bin/env python3
import subprocess as sp
import sys
import os
import warnings


rule pycoQc_fastq:
    input:
        rename = "{sample_id}_renamed.done"
    output:
        fastqc = temp("{sample_id}_qc.done")
    run:
        this_sample = config["data"][wildcards.sample_id]
        if not os.path.exists("FASTQC_Project_"+this_sample["Sample_Project"]):
            os.mkdir("FASTQC_Project_"+this_sample["Sample_Project"])
        if this_sample["barcode_kits"] is not "no_bc":
            bs = this_sample["index_id"]
            bs = [s for s in " ".join(list(bs)).split() if s.isdigit()]
            bs = "".join(bs)
            print("bs:", bs)
            barcode = "barcode"+bs
            fastqc = os.path.join(os.path.join("FASTQC_Project_"+this_sample["Sample_Project"],
                                               "Sample_"+this_sample["Sample_ID"]))
            if not os.path.exists(fastqc):
                os.mkdir(fastqc)

                cmd = config["pycoQc"]["pycoQc"]+" "
                cmd += "fastq/sequencing_summary_"+barcode+".txt "
                cmd += "-o "+fastqc+"/pycoqc_"+barcode+".html"
                print(cmd)
                sp.check_call(cmd, shell=True)
            else:
                warnings.warn("This directory: ",fastqc, "already exists, computing qc is skipped "
                              "and the folder stays untouched.")

            sp.check_call("touch "+output.fastqc, shell=True)
        else:
            assert(len(config["data"])==1)
            project_name = list(config["data"].values())[0]["Sample_Project"]
            fastqc = os.path.join("FASTQC_Project_"+project_name)
            print(fastqc)
            assert(os.path.exists("FASTQC_Project_"+project_name) == 0)
            os.mkdir("FASTQC_Project_"+project_name)
            cmd = config["pycoQc"]["pycoQc"]+" "
            cmd += "fastq/sequencing_summary.txt "
            cmd += "-o "+fastqc+"/pycoqc.html"
            print(cmd)
            sp.check_call(cmd, shell=True)


rule pycoQc_bam:
    input:

    output:

    run:
        if bc_kit  == 'no_bc':
            assert(len(data.items()) == 1)
            sample_project = list(data.values())[0]["Sample_Project"]
            sample_id = list(data.values())[0]["Sample_ID"]
            bam = list(data.values())[0]["Sample_Name"]+".bam"
            group=sample_project.split("_")[-1].lower()
            final_path = config["paths"]["groupDir"]+group+"/sequencing_data/OxfordNanopore/"+config["input"]["name"]
            path_to_bam = final_path+"/Analysis_"+sample_project+"/mapping_on_"+ref+"/"+sample_id
            bam = os.path.join(path_to_bam, bam)
            cmd = config["pycoQc"]["pycoQc"]+" "
            cmd += wd+"/fastq/sequencing_summary.txt "
            cmd += " -a "+bam
            cmd += " -o "+path_to_bam+"/pycoqc_bam.html"
            print(cmd)
            sp.check_call(cmd, shell=True)
#     else:
#         for k, v in data.items():
#             group=v["Sample_Project"].split("_")[-1].lower()
#             final_path = config["paths"]["groupDir"]+group+"/sequencing_data/OxfordNanopore/"+config["input"]["name"]
#             path_to_bam = final_path+"/Analysis_"+v["Sample_Project"]+"/mapping_on_"+ref+"/"+v["Sample_ID"]
#             bam = v["Sample_Name"]+".bam"
#             bam = os.path.join(path_to_bam, bam)
#             barcode = v["index_id"]
#             cmd = config["pycoQc"]["pycoQc"]+" "
#             cmd += wd+"/sequencing_summary_"+barcode+".txt "
#             cmd += " -a "+bam
#             cmd += " -o "+path_to_bam+"/pycoqc_"+v["Sample_Name"]+"_"+barcode+".html"
#             print(cmd)
#             sp.check_call(cmd, shell=True)
