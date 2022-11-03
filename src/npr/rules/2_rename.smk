#!/usr/bin/env python3
import os
import sys
import subprocess as sp
import warnings
from pathlib import Path
from npr.snakehelper import config_to_splitseqsummary
import shutil
import glob
import pandas as pd

# create a pandas dataframe of samples to get the sample : project relationship
metadata = dict(config["data"])
del metadata["projects"]
del metadata["samples"]
metadata = pd.DataFrame(metadata).T

rule rename:
    input:
        "flags/1_basecalling.done"
    output:
        flag=touch("flags/2_renamed.done"),
        fastqs=expand_project_path("Project_{project}/Sample_{sample_id}/pass/{sample_name}.fastq.gz"),  
        summaries=expand_project_path("Project_{project}/Sample_{sample_id}/sequencing_summary.txt"),
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
                pass_out = os.path.join(
                    sampleid_dir,
                    'pass',
                    samDic['Sample_Name'] + '.fastq.gz'
                )
                fail_out = os.path.join(
                    sampleid_dir,
                    'fail',
                    samDic['Sample_Name'] + '.fastq.gz'
                )
                logfile.write("Creating directories\n")
                os.mkdir(project_dir)
                os.mkdir(sampleid_dir)
                os.mkdir(
                    os.path.join(
                        sampleid_dir, 'pass'
                    )
                )
                os.mkdir(
                    os.path.join(
                        sampleid_dir, 'fail'
                    )
                )
                logfile.write("Passing fastq files.\n")
                passlist = glob.glob(
                    os.path.join('fastq','pass','*fastq.gz')
                )
                cmd = ['cat'] + passlist
                with open(pass_out, 'w') as f:
                    sp.call(cmd, stdout=f)
                for f in passlist:
                    os.remove(f)
                #fail
                logfile.write("failing fastq files.\n")
                faillist = glob.glob('fastq/fail/*fastq.gz')
                cmd = ['cat'] + faillist
                with open(fail_out, 'w') as f:
                    sp.call(cmd, stdout=f)
                for f in faillist:
                    os.remove(f)
                logfile.write("Copy sequencing_summary to sample folder.\n")
                shutil.copy(
                    'fastq/sequencing_summary.txt',
                    os.path.join(
                        sampleid_dir, 
                        'sequencing_summary.txt'
                )
            else:
                # To Do - make pass/fail splits here as well.
                logfile.write("barcoding detected.\n")
                splitcmd = config_to_splitseqsummary(config)
                logfile.write("splitcmd: {}\n".format(splitcmd))
                sp.check_call(splitcmd, shell=True)
                for sample_id in config['data']['samples']:
                    logfile.write("Renaming sample {}\n".format(sample_id))
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
                    # Cat only the specific sample ofcourse..
                    cmd = ['cat']
                    for fqFile in glob.glob('fastq/pass/{}/*fastq.gz'.format(samDic['index_id'])):
                        cmd.append(fqFile)
                    for fqFile in glob.glob('fastq/fail/{}/*fastq.gz'.format(samDic['index_id'])):
                        cmd.append(fqFile)
                    logfile.write("Running fastq cat for barcode {}.\n".format(samDic['index_id']))
                    with open(fq_out, 'w') as f:
                        sp.call(cmd, stdout=f)
                    shutil.copy(
                        'sequencing_summary_{}.txt'.format(samDic['index_id']),
                        os.path.join(
                            sampleid_dir,
                            'sequencing_summary_{}.txt'.format(samDic['Sample_Name'])
                        )
                    )
                shutil.move(
                    'sequencing_summary_unclassified.txt',
                    'reports/sequencing_summary_unclassified.txt'
                )
                
