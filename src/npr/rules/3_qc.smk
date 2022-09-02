#!/usr/bin/env python3
import subprocess as sp
import sys
import os
import warnings
import glob
from npr.snakehelper import config_to_pycoqc
from npr.snakehelper import grab_seqsummary


rule pycoQc_fastq:
    input:
        "flags/2_renamed.done"
    output:
        touch('flags/3_qc.done')
    log:
        log = 'log/3_qc.log'
    run:
        # Assumption is only 1 project per flowcell..
        with open(log.log, 'w') as logfile:
            for sample_id in config['data']['samples']:
                logfile.write("pycoqc on sample {}\n".format(sample_id))
                sample_dic = config['data'][sample_id]
                sample_project = sample_dic['Sample_Project']
                fqc_dir = os.path.join(
                    'FASTQC_Project_' + sample_project
                )
                fqc_sampledir = os.path.join(
                    fqc_dir,
                    'Sample_' + sample_id
                )
                project_sampledir = fqc_sampledir.replace('FASTQC_', '')
                if not os.path.exists(fqc_dir):
                    os.mkdir(fqc_dir)
                os.mkdir(fqc_sampledir)
                cmd = config_to_pycoqc(
                    config,
                    grab_seqsummary(project_sampledir),
                    fqc_sampledir,
                    sample_id,
                    config['bc_kit']
                )
                print(cmd)
                logfile.write("pycoqc command:\n")
                logfile.write("{}\n".format(cmd))
                sp.check_call(cmd, shell=True)