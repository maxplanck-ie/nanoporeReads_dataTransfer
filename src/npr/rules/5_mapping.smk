#!/usr/bin/env python3
import subprocess as sp
import os
import glob
import shutil
from npr.snakehelper import get_seqdir
from npr.snakehelper import config_to_mapcmd
from npr.snakehelper import grab_seqsummary


rule dna_mapping:
    input:
        'flags/4_transfer.done'
    output:
        touch('flags/5_mapping.done')
    log:
        log = 'log/5_mapping.log'
    run:
        with open(log.log, 'w') as logfile:
            # Get final directory.
            p = config['data']['projects'][0]
            project_dir = 'Project_' + p
            fqc_dir = 'FASTQC_' + project_dir
            group = project_dir.split("_")[-1].lower()
            logfile.write("Transfer: {} & {}\n".format(
                project_dir,
                fqc_dir
            ))
            groupdir = os.path.join(
                config["paths"]["groupDir"],
                group
            )
            if not os.path.exists(groupdir): # If external (no /data/pi/ path exists)
                group_path = os.path.join(
                    config["paths"]["external_groupDir"],
                    group,
                    "sequencing_data/OxfordNanopore"
                )
                if not os.path.exists(group_path):
                    os.makedirs(group_path)
            else:
                group_path = get_seqdir(
                    groupdir, "sequencing_data"
                )
            final_path = os.path.join(group_path, config["input"]["name"])
            if config['info_dict']['organism'] in list(config['genome'].keys()):
                genome = config['info_dict']['organism']
                logfile.write("Valid organism: {}\n".format(genome))
                pref, post = config_to_mapcmd(config)
                if pref:
                    analysis_dir = os.path.join(
                        final_path,
                        'Analysis_' + p
                    )
                    if os.path.exists(analysis_dir):
                        shutil.rmtree(analysis_dir)
                    logfile.write("Making {}\n".format(analysis_dir))
                    os.mkdir(analysis_dir)
                    genome_analysis_dir = os.path.join(
                        analysis_dir,
                        config['info_dict']['protocol'] + '_' + genome
                        )
                    os.mkdir(genome_analysis_dir)
                    for sample in config['data']['samples']:
                        pref, post = config_to_mapcmd(config)
                        logfile.write("Sample {}\n".format(sample))
                        logfile.write("Sample in config: {}\n".format(
                            config['data'][sample]['Sample_Name']
                        ))
                        fqFile = os.path.join(
                            final_path,
                            'Project_' + p,
                            'Sample_' + sample,
                            config['data'][sample]['Sample_Name'] + '.fastq.gz'
                        )
                        logfile.write("fqFile path = {}".format(fqFile))
                        pref.append(fqFile)
                        pref.append('|')
                        bamfile = os.path.join(
                                genome_analysis_dir,
                                config['data'][sample]['Sample_Name'] + '.bam'
                            )
                        post.append(
                            bamfile
                        )
                        #postsample.append('-')
                        mapCmd = pref+post
                        logfile.write(
                            "Mapping command invoked:\n"
                        )
                        logfile.write(
                            ' '.join(mapCmd) + "\n"
                        )
                        # Run alignment.
                        sp.check_call(' '.join(mapCmd), shell=True)
                        # index bam file.
                        ixCmd = [
                            'samtools',
                            'index',
                            '-@',
                            '25',
                            bamfile
                        ]
                        sp.check_call(' '.join(ixCmd), shell=True)
                        # pycoQc
                        qcCmd = [
                            'pycoQC',
                            '--summary_file',
                            grab_seqsummary(
                                os.path.join(
                                    final_path,
                                    project_dir,
                                    'Sample_' + sample
                                )
                            ),
                            '-a',
                            bamfile,
                            '-o',
                            bamfile.replace('.bam', '.html')
                        ]
                        sp.check_call(' '.join(qcCmd), shell=True)
