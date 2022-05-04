#snakefile
import os
import yaml
import pandas as pd
import warnings
import subprocess as sp
import misc.sendEmail as email

# f = open("config.yaml")
# globals().update(yaml.load(f))
# f.close()

# include paths:
include: os.path.join("rules", "baseCalling.Snakefile")
include: os.path.join("rules", "qc.py")
include: os.path.join("rules", "rename.Snakefile")
include: os.path.join("rules", "data_transfer.py")
include: os.path.join("rules", "mapping.py")


def run_basecalling():
    return [expand("demux.done")]


def run_bc_split():
    return [expand("bc.split")]


def run_rename():
    return ["renamed.done"]


def run_pycoqc():
    return ["fastqQC.done"]


def run_data_transfer():
    return ["transfer.done"]


def run_mapping():
    return ["mapping.done"]


def run_bamqc():
    return ["bamQC.done"]


def qc2deepseq():
    return


def transfer2rapidus():
    return

# def run_contamination_report():
#     report_contamination(config, data, args.protocol)


rule all:
    input:
        run_basecalling(),
        run_bc_split(),
        run_rename(),
        run_pycoqc(),
        run_data_transfer(),
        run_mapping(),
        run_bamqc()
        # # run_bamCompare()
        # run_contamination_report()


onsuccess:
    shell("touch analysis.done")
onerror:
    email.sendEmail("snakemake pipeline failed!", "report a crash", config['email']['from'], config['email']['to'],
                    config['email']['host'])
