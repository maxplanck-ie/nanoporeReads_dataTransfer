#!/usr/bin/env python3

"""
Warning: 
This workflow assumes that all barcodes are present and subsequent fastqs can be generated.
This is particularly problematic for missing barcodes.
Currently the rule will fail after waiting for missing file.
Do not touch empty files as subsequent QC step will fail with empty input.
Instead the SampleSheet (config['data']) should be adjusted
"""

import os
import sys
import subprocess as sp
import warnings
from pathlib import Path
from npr.snakehelper import config_to_splitseqsummary
import shutil
import glob
import pandas as pd
import traceback

# create a pandas dataframe of samples to get the sample : project relationship
metadata = dict(config["data"])
del metadata["projects"]
del metadata["samples"]
metadata = pd.DataFrame(metadata).T

def cat_list(file_list, outfile, msg=""):
    """
    concatenate files (fastq.gz) in file_list and send them to outfile
    include checks in case list is empty
    """

    # file_list might be empty (e.g if barcode is not found)
    if file_list:
        cmd = ['cat'] + file_list
        with open(outfile, 'w') as f:
            sp.call(cmd, stdout=f)
    else:
        print('file_list is empty. {}'.format(msg))
        sys.exit(1)
        #Path(outfile).touch()


rule rename_files:
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
                if not os.path.exists(project_dir):
                    os.mkdir(project_dir)
                if not os.path.exists(sampleid_dir):
                    os.mkdir(sampleid_dir)
                logfile.write("Creating pass.")
                passdir = os.path.join(
                    sampleid_dir, 'pass'
                )
                faildir = os.path.join(
                    sampleid_dir, 'fail'
                )
                if not os.path.exists(passdir):
                    os.mkdir(passdir)
                if not os.path.exists(faildir):
                    os.mkdir(faildir)
                logfile.write("Passing fastq files.\n")
                passlist = glob.glob(
                    os.path.join('fastq','pass','*fastq.gz')
                )
                cat_list(passlist, pass_out, msg="No fail/*fastq.gz")
                #for f in passlist:
                #    os.remove(f)
                #fail
                logfile.write("failing fastq files.\n")
                faillist = glob.glob('fastq/fail/*fastq.gz')
                cat_list(faillist, fail_out, msg="No fail/*fastq.gz")
                #for f in faillist:
                #    os.remove(f)

                logfile.write("Copy sequencing_summary to sample folder.\n")
                ss_file = 'fastq/sequencing_summary.txt'
                ss_new = os.path.join(sampleid_dir, 'sequencing_summary.txt')
                if Path(ss_file).exists():
                    shutil.copy( ss_file, ss_new )
                    logfile.write('copied {} to {}'.format(ss_file, ss_new))
                else:
                    #Path(ss_new).touch()
                    logfile.write('no such file {}. Check basecalling.'.format(ss_file))
            else:
                logfile.write("barcoding detected.\n")
                splitcmd = config_to_splitseqsummary(config)
                logfile.write("splitcmd: {}\n".format(splitcmd))
                # Split fastq/sequencing_summary.txt by barcode
                sp.check_call(splitcmd, shell=True)
                for sample_id in config['data']['samples']:
                    logfile.write("Renaming sample {}\n".format(sample_id))
                    samDic = config['data'][sample_id]
                    # Some exceptions...
                    samDic['index_id'] = samDic['index_id'].replace('021', '21').replace('022', '22').replace('023', '23').replace('024', '24')
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
                    logfile.write("Creating directories for {}\n".format(sample_id))
                    if not os.path.exists(project_dir):
                        os.mkdir(project_dir)
                    if not os.path.exists(sampleid_dir):
                        os.mkdir(sampleid_dir)
                    logfile.write("Creating pass.")
                    passdir = os.path.join(
                        sampleid_dir, 'pass'
                    )
                    faildir = os.path.join(
                        sampleid_dir, 'fail'
                    )
                    if not os.path.exists(passdir):
                        os.mkdir(passdir)
                    logfile.write("Passing fastq files.\n")
                    passlist = glob.glob(
                        os.path.join('fastq','pass','{}'.format(samDic['index_id']), '*fastq.gz')
                    )

                    cat_list(passlist, pass_out, msg="No pass/" + samDic['index_id'])


                    #for f in passlist:
                    #    os.remove(f)
                    #fail
                    logfile.write("failing fastq files.\n")
                    faillist = glob.glob(
                        os.path.join('fastq','fail','{}'.format(samDic['index_id']), '*fastq.gz')
                    )
                    if not os.path.exists(faildir):
                        os.mkdir(faildir)

                    cat_list(faillist, fail_out, msg="No fail/" + samDic['index_id'])

                    ss_file = 'sequencing_summary_{}.txt'.format(samDic['index_id'])
                    ss_new = os.path.join(sampleid_dir, 'sequencing_summary.txt')
                    if Path(ss_file).exists():
                        shutil.copy( ss_file, ss_new )
                        logfile.write('copied {} to {}'.format(ss_file, ss_new))
                    else:
                        #Path(ss_new).touch()
                        logfile.write('no such file {}'.format(ss_file))
                        logfile.write('Some barcodes from sample sheet seem to be missing in data: {}'.format(ss_file))
