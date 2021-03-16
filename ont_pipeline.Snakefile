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


def run_basecalling():
    return [expand("fastq")]


def run_pycoqc():
    file_list = []
    file_list.append([expand("{sample_id}_qc.done", sample_id = config["data"].keys())])
    return file_list


def run_data_transfer():
    file_list = []
    file_list.append([expand("{sample_id}.transferred", sample_id = config["data"].keys())])
    return file_list


# def run_mapping():
#     if ("RNA" in config["info_dict"]["kit"]) or (args.protocol == 'rna') or (args.protocol == 'cdna'):
#         print("minimap2 starts")
#         mapping_rna(config, data, args.reference)
#         pycoQc_bam(config,data, bc_kit, args.reference)
#     else:
#         mapping_dna(config, data, args.reference)
#         pycoQc_bam(config,data, bc_kit, args.reference)
#     print("data has been mapped")
#
#
# def run_contamination_report():
#     report_contamination(config, data, args.protocol)


rule all:
    input:
        run_basecalling(),
        expand("bc.split"),
        run_pycoqc()
        run_data_transfer()
        # run_mapping()
        # run_contamination_report()



onsuccess:
    shell("rm bs.split")
# onerror:
