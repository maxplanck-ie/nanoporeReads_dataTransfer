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
        fastq="Project_{project}/Sample_{sample_id}/pass/{sample_name}_porechop.fastq.gz",
        info="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_porechop.info",
    threads: 8
    params:
        # guppy returns "U" for RNA so -abi will not work
        flag="-abi" if config['info_dict']['protocol'] != 'rna' else "",
    log:
        err='log/porechop/project-{project}_id-{sample_id}_name-{sample_name}.err'
    shell:'''
        porechop_abi {params.flag} -t {threads} -i {input.fastq} -o {output.fastq} > {output.info} 2> {log.err}
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

rule qc_kraken:
    input:
        fastq="Project_{project}/Sample_{sample_id}/pass/{sample_name}.fastq.gz"
    output:
        report="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_kraken.report",
        output="FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_kraken.output",
    params:
        db=config['kraken']['db'],
        threads=config['kraken']['threads'],
    log:
        stderr="log/kraken/project-{project}_id-{sample_id}_name-{sample_name}.log"
    shell:'''
        kraken2 -db {params.db} --threads {params.threads} --use-names {input.fastq}  --output {output.output} --report {output.report} 2> {log.stderr}
    '''

rule qc_multiqc: 
    input: 
        pycoQc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_pycoqc.html"),
        porechopQc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_porechop.info"),
        fastqc=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_fastqc.html"),
        krakenQC=expand_project_path("FASTQC_Project_{project}/Sample_{sample_id}/{sample_name}_kraken.report"),
    output: 
        html="FASTQC_Project_{project}/multiqc/multiqc_report.html", 
    params:
        pdir="FASTQC_Project_{project}",
        mdir="FASTQC_Project_{project}/multiqc",
    log:
        err="log/multiqc/project-{project}.err"
    shell:'''
        multiqc -o {params.mdir} {params.pdir} 2> {log.err}
    '''

rule qc_finalize: 
    input: 
        multiqc=expand_project_path("FASTQC_Project_{project}/multiqc/multiqc_report.html"),
    output: 
        touch("flags/3_qc.done"),