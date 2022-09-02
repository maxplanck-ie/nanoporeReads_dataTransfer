#!/usr/bin/env python3

import subprocess as sp
import sys
import os
import warnings
import glob
import shutil
from npr.snakehelper import get_seqdir
from npr.snakehelper import genome_index
from npr.snakehelper import overwrite_dir

rule data_transfer:
    input:
        'flags/3_qc.done'
    output:
        touch('flags/4_transfer.done')
    log:
        log = 'log/4_transfer.log'
    run:
        with open(log.log, 'w') as logfile:
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
            logfile.write("group_path = {}\n".format(group_path))
            final_path = os.path.join(group_path, config["input"]["name"])
            config['data']['finalpath'] = final_path
            if not os.path.exists(final_path):
                os.mkdir(final_path)
            # copy over data.
            logfile.write("init writing.\n")
            cp_proj = overwrite_dir(project_dir, final_path)
            logfile.write("project = {}\n".format(cp_proj))
            cp_pod5 = overwrite_dir('pod5', os.path.join(final_path, project_dir))
            logfile.write("pod5 = {}\n".format(cp_pod5))
            cp_fqc = overwrite_dir(fqc_dir, final_path)
            logfile.write("pod5 = {}\n".format(cp_fqc))