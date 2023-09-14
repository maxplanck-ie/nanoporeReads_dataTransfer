#!/usr/bin/env python3
import subprocess as sp
import sys
import os
import warnings
import glob
from npr.snakehelper import config_to_pycoqc
from npr.snakehelper import grab_seqsummary

rule qc_porechop:
    input:
        fastq="Project_{project}/Sample_{sample_id}/pass/{sample_name}.fastq.gz",
    output:
        fastq="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_porechop.fastq",
        info="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_porechop.info",
        summary="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_porechop.summary",
    log:
        err='log/porechop/project-{project}_id-{sample_id}_name-{sample_name}.err'
    shell:'''
        porechop_abi -t {threads} -i {input.fastq} -o {output.fastq} > {output.info} 2> {log.err}
        # || avoids non-zero exit status if grep returns nothing
        # see: https://unix.stackexchange.com/questions/330660
        ( grep "adapter" {output.info} || [[ $? == 1 ]] )  > {output.summary}
        sed -n '/Barcode  Reads/,$p'  {output.info} >> {output.summary} 2> {log.err}
        '''

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
    params:
        odir = lambda wildcards: "FASTQC_Project_{}/Sample_{}/".format(wildcards.project, wildcards.sample_id)
    threads: 16
    log:
        "log/fastqc/project-{project}_id-{sample_id}_name-{sample_name}.log"
    shell:'''
    fastqc --memory=4096 -t {threads} -o {params.odir} --dir {params.odir} --quiet {input}
    '''

rule qc_reformat: 
    input: 
        pycoQc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_pycoqc.html"),
        porechopQc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_porechop.summary"),
        fastqc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_fastqc.html"),
    output: 
        touch("flags/3_qc.done")
