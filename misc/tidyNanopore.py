#!/usr/bin/env python

# Tidy up Nanopore directories, it will remove subdirectories for flow cells that have been archived

import os
import glob
import shutil
import argparse
import subprocess
import logging

parser = argparse.ArgumentParser(description='Tidy Nanopore directories that have been archived.')
parser.add_argument('-d', '--delete',
                    action='store_true',
                    help='Delete the directories (otherwise, just show the ones to be removed)')
parser.add_argument('-t', '--target-dir',
                    help='Target directory',
                    required=True)
parser.add_argument('-l', '--last-run',
                    help='Date of last run to be considered, format YYYYMMDD',
                    required=True)
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    help='Verbose mode')

args = parser.parse_args()

do_delete  = args.delete
target_dir = args.target_dir
last_run   = args.last_run

if args.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

# Get all of the flow cells
dirs = glob.glob(target_dir + "/*")

for d in sorted(dirs):
    if not os.path.isdir(d):
        continue

    # Skip directories that have NOT been archived
    bak_dir = d + ".is_on_archive"
    if not os.path.exists(bak_dir):
        logging.debug(f"{d} is not on archive, skipping")
        continue

    # Skip younger directories
    date_dir = os.path.basename(d)[:8]
    if date_dir > last_run:
        logging.debug(f"{d} is younger than {last_run}, skipping")
        continue

    # Check if the flowcell has been analyzed
    if not os.path.exists(d + "/analysis.done"):
        logging.debug(f"{d} has not been analyzed, skipping")
        continue

    # define large subdirectories to be removed
    subdirs = glob.glob("{}/Project_*".format(d)) # legacy
    subdirs.extend(glob.glob("{}/fast5*".format(d)))
    subdirs.extend(glob.glob("{}/fastq*".format(d))) # legacy
    subdirs.extend(glob.glob("{}/pod5*".format(d)))
    subdirs.extend(glob.glob("{}/bam*".format(d)))
    subdirs.extend(glob.glob("{}/transfer/Project_*/Data".format(d)))
    subdirs.extend(glob.glob("{}/transfer/Project_*/Analysis*".format(d)))

    # check if subdir is actually a directory
    for d2 in subdirs:
        if not os.path.isdir(d2):
            subdirs.remove(d2)

    # define txt files to be compressed (legacy)
    compress = glob.glob("{}/*.txt".format(d))
    # check if txt file is actually a file
    for f in compress:
        if not os.path.isfile(f):
            compress.remove(f)

    # do nothing if no subdirs
    if len(subdirs) + len(compress) == 0:
        logging.debug(f"{d} has no subdirs to be deleted or removed, skipping")
        continue

    if do_delete:
        logging.info('Tidy ' + d)
        for d2 in subdirs:
            logging.info("Removing {}".format(d2))
            shutil.rmtree("{}".format(d2))
        for f in compress:
            logging.info("Compressing {}".format(f))
            subprocess.check_call(['gzip', f])
    else:
        for d2 in subdirs:
            logging.info("Would removed: {}".format(d2))
        for f in compress:
            logging.info("Would compress {}".format(f))

logging.info("Cleaning done.")
