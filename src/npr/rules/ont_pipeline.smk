#snakefile
import os
import yaml
import pandas as pd
import warnings
import subprocess as sp
from npr.communication import send_email

include: os.path.join(config['paths']['rulesPath'], 'basecalling.smk')
include: os.path.join(config['paths']['rulesPath'], "qc.py")
include: os.path.join(config['paths']['rulesPath'], "rename.smk")
include: os.path.join(config['paths']['rulesPath'], "data_transfer.py")
include: os.path.join(config['paths']['rulesPath'], "mapping.py")


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

rule all:
    input:
        run_basecalling(),
        run_bc_split(),
        run_rename(),
        run_pycoqc(),
        run_data_transfer(),
        run_mapping(),
        run_bamqc()
