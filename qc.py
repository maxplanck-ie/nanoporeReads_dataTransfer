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

def pycoQc_fastq(config, data, bc_kit):
    wd = os.path.join(config["info_dict"]["flowcell_path"])
    if bc_kit is not "no_bc":
        cmd = config["pycoQc"]["barcodeSplit"]
        cmd += wd+"/fastq/sequencing_summary.txt "
        cmd += "-o "+wd
        print(cmd)
        sp.check_call(cmd, shell=True)
        for k, v in data.items():
            fastqc = os.path.join(wd+"/FASTQC_Project_"+v["Sample_Project"])
            if not os.path.exists(wd+"/FASTQC_Project_"+v["Sample_Project"]):
                os.mkdir(wd+"/FASTQC_Project_"+v["Sample_Project"])
            barcode = v["index_id"]
            os.mkdir(fastqc+"/Sample_"+v["Sample_ID"])
            fastqc = os.path.join(fastqc+"/Sample_"+v["Sample_ID"])
            cmd = config["pycoQc"]["pycoQc"]+" "
            cmd += wd+"/sequencing_summary_"+barcode+".txt "
            cmd += "-o "+fastqc+"/pycoqc_"+barcode+".html"
            print(cmd)
            sp.check_call(cmd, shell=True)

        cmd = config["pycoQc"]["pycoQc"]+" "
        cmd += wd+"/fastq/sequencing_summary.txt "
        cmd += "-o "+wd+"/pycoqc.html"
        print(cmd)
        sp.check_call(cmd, shell=True)
    else:
        assert(len(data.items()) == 1)
        project_name = list(data.values())[0]["Sample_Project"]
        fastqc = os.path.join(wd+"/FASTQC_Project_"+project_name)
        print(fastqc)
        assert(os.path.exists(wd+"/FASTQC_Project_"+project_name) == 0)
        os.mkdir(wd+"/FASTQC_Project_"+project_name)
        cmd = config["pycoQc"]["pycoQc"]+" "
        cmd += wd+"/fastq/sequencing_summary.txt "
        cmd += "-o "+fastqc+"/pycoqc.html"
        print(cmd)
        sp.check_call(cmd, shell=True)


def pycoQc_bam(config, data, bc_kit, ref):
    wd = os.path.join(config["info_dict"]["flowcell_path"])
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
