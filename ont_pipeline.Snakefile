#snakefile
import os
import yaml
import pandas as pd

# f = open("config.yaml")
# globals().update(yaml.load(f))
# f.close()

# include paths:
include: os.path.join("rules", "baseCalling.py")
include: os.path.join("rules", "qc.py")
include: os.path.join("rules", "rename.py")
include: os.path.join("rules", "data_transfer.py")
include: os.path.join("rules", "mapping.py")


def run_basecalling():
    return "fastq"


def run_pycoqc():
    file_list = []
    file_list.append([expand("{sample_id}_qc.done", sample_id = config["data"].keys())])
    return file_list


def run_data_transfer():
    file_list = []
    file_list.append([expand("{sample_id}.transferred", sample_id = config["data"].keys())])
    return file_list


def run_mapping():
    file_list = []
    file_list.append([expand("{sample_id}.bam", sample_id = config["data"].keys())])
    return file_list


# def run_contamination_report():
#     report_contamination(config, data, args.protocol)


rule all:
    input:
        run_basecalling(),
        expand("bc.split"),
        run_pycoqc(),
        run_data_transfer(),
        run_mapping(),
        # run_bam_qc(),
        # run_contamination_report()



onsuccess:
   shell("rm bs.split")
# onerror:
