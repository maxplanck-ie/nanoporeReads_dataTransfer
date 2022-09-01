#!/usr/bin/env python3
import os
import sys
import subprocess as sp
import warnings
from pathlib import Path
from npr.snakehelper import config_to_splitseqsummary
import shutil


rule rename:
    input:
        "flags/1_basecalling.done"
    output:
        touch("flags/2_renamed.done")
    log:
        log = 'log/2_rename.log'
    run:
        with open(log.log, 'w') as logfile:
            if config['bc_kit'] == 'no_bc':
                logfile.write("no barcoding, 1 sample only?\n")
                # Make sure there is only 1 sample.
                assert(len(config['data']['samples']) == 1)
                logfile.write("1 sample only asserted..\n")
                logfile.write("Only 1 sample.\n")
                sample_id = config['data']['samples'][0]
                samDic = config['data'][sample_id]
                project_dir = "Project_" + samDic['Sample_Project']
                sampleid_dir = os.path.join(
                    project_dir,
                    'Sample_' + sample_id
                )
                fq_out = os.path.join(
                    sampleid_dir,
                    samDic['Sample_Name'] + '.fastq.gz'
                )
                logfile.write("Creating directories\n")
                os.mkdir(project_dir)
                os.mkdir(sampleid_dir)
                cmd = [
                    'cat'
                ]
                for fqFile in glob.glob('fastq/pass/*fastq.gz'):
                    cmd.append(fqFile)
                for fqFile in glob.glob('fastq/fail/*fastq.gz'):
                    cmd.append(fqFile)
                logfile.write("Running fastq cat.\n")
                with open(fq_out, 'w') as f:
                    sp.call(cmd, stdout=f)
                logfile.write("Copy sequencing_summary to sample folder.\n")
                shutil.copy(
                    'fastq/sequencing_summary.txt',
                    os.path.join(sampleid_dir, 'sequencing_summary.txt')
                )
            else:
                logfile.write("barcoding detected.\n")
                splitcmd = config_to_splitseqsummary(config)
                sp.check_call(splitcmd)
                for sample_id in config['data']:
                    print("rename some barcodes & crap")