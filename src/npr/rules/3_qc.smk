#!/usr/bin/env python3
import subprocess as sp
import sys
import os
import warnings
import glob
from npr.snakehelper import config_to_pycoqc
from npr.snakehelper import grab_seqsummary


rule qc_pycoqc:
    input:
        fastq="Project_{project}/Sample_{sample_id}/pass/{sample_name}.fastq.gz",
        sequencing_summary="Project_{project}/Sample_{sample_id}/sequencing_summary.txt",
    output:
        "FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_pycoqc.html"
    log:
        'log/pycoqc/project-{project}_id-{sample_id}_name-{sample_name}.log'
    shell:
        "pycoQC --summary_file {input.sequencing_summary} " 
        "{config[pycoQc][pycoQc_opts]} "
        "-o {output} > {log} 2>&1"

rule qc_fastqc: 
    input: 
        "Project_{project}/Sample_{sample_id}/pass/{sample_name}.fastq.gz"
    output: 
        html="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_fastqc.html",
        zip="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_fastqc.zip", # suffix _fastqc.zip is necessary for multiqc
    params: "--quiet"
    threads: 16
    log:
        "log/fastqc/project-{project}_id-{sample_id}_name-{sample_name}.log"
    wrapper: 
        "v1.12.2/bio/fastqc"
    

rule qc_reformat: 
    input: 
        pycoQc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_pycoqc.html"),
        fastqc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_fastqc.html"),
    output: 
        touch("flags/3_qc.done")
