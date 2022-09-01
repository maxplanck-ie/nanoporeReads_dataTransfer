#!/usr/bin/env python3
import subprocess as sp
import os
import glob
import shutil
from npr.snakehelper import get_seqdir
from npr.snakehelper import config_to_mapcmd


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
                        'Analysis_' + p + '_' + '-'.join(
                            [
                                genome,
                                config['info_dict']['protocol']
                            ]
                        )
                    )
                    if os.path.exists(analysis_dir):
                        shutil.rmtree(analysis_dir)
                    logfile.write("Making {}\n".format(analysis_dir))
                    os.mkdir(analysis_dir)
                    os.mkdir(
                        os.path.join(
                            analysis_dir,
                            'bamfiles'
                        )
                    )
                    for sample in config['data']['samples']:
                        prefsample = pref
                        postsample = post
                        fqFile = os.path.join(
                            final_path,
                            'Project_' + p,
                            'Sample_' + sample,
                            config['data'][sample]['Sample_Name'] + '.fastq.gz'
                        )
                        prefsample.append(fqFile)
                        prefsample.append('|')
                        postsample.append(
                            os.path.join(
                                analysis_dir,
                                'bamfiles',
                                config['data'][sample]['Sample_Name'] + '.bam'
                            )
                        )
                        mapCmd = prefsample+postsample
                        logfile.write(
                            "Mapping command invoked:\n"
                        )
                        logfile.write(
                            ' '.join(mapCmd) + "\n"
                        )
                        print(' '.join(mapCmd))
                        sp.check_call(mapCmd)
