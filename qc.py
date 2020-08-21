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

def pycoQc(config, data):
    wd = os.path.join(config["info_dict"]["flowcell_path"])
    cmd = config["pycoQc"]["barcodeSplit"]
    cmd += wd+"/fastq/sequencing_summary.txt "
    cmd += "-o "+wd
    print(cmd)
    QC_path = ""
    sp.check_call(cmd, shell=True)
    for k, v in data.items():
        fastqc = os.path.join(wd+"/FASTQC_Project_"+v["Sample_Project"])
        if not os.path.exists(wd+"/FASTQC_Project_"+v["Sample_Project"]):
           os.mkdir(wd+"/FASTQC_Project_"+v["Sample_Project"])
        QC_path = fastqc
        if v["barcode_kits"] is not "no_bc":
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
    cmd += "-o "+QC_path+"/pycoqc.html"
    print(cmd)
    sp.check_call(cmd, shell=True)
