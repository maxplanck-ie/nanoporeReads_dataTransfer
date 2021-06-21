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
    return [expand("demux.done")]


def run_bc_split():
    return [expand("bc.split")]

def run_rename():
    file_list = []
    file_list.append([expand("{sample_id}_renamed.done", sample_id = config["data"].keys())])
    return file_list


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
    file_list.append([expand("{sample_id}.mapped", sample_id = config["data"].keys())])
    return file_list


def run_bamqc():
    file_list = []
    file_list.append([expand("{sample_id}.bam.qc", sample_id = config["data"].keys())])
    return file_list


# def run_contamination_report():
#     report_contamination(config, data, args.protocol)


rule all:
    input:
        "demux.done",
        "bc.split",
        run_rename(),
        run_pycoqc(),
        run_data_transfer(),
        run_mapping(),
        run_bamqc(),
        # run_bamCompare()
        # run_contamination_report()



# onsuccess:
#     shell("rm demux.done bc.split *.transferred *renamed.done *qc.done *.mapped *bam.qc")
# onerror:
