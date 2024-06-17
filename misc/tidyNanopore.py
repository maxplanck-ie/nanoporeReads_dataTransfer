#!/usr/bin/env python

# Tidy up Nanopore directories, it will remove subdirectories for flow cells that have been archived

import argparse
import glob
import logging
import os
import shutil
import subprocess

parser = argparse.ArgumentParser(
    description="Tidy Nanopore directories that have been archived."
)
parser.add_argument(
    "-d",
    "--delete",
    action="store_true",
    help="Delete the directories (otherwise, just show the ones to be removed)",
)
parser.add_argument("-t", "--target-dir", help="Target directory", required=True)
parser.add_argument(
    "-l",
    "--last-run",
    help="Date of last run to be considered, format YYYYMMDD",
    required=True,
)
parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")

args = parser.parse_args()

do_delete = args.delete
target_dir = args.target_dir
last_run = args.last_run

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
    subdirs = glob.glob(f"{d}/Project_*")  # legacy
    subdirs.extend(glob.glob(f"{d}/fast5*"))
    subdirs.extend(glob.glob(f"{d}/fastq*"))  # legacy
    subdirs.extend(glob.glob(f"{d}/pod5*"))
    subdirs.extend(glob.glob(f"{d}/bam*"))
    subdirs.extend(glob.glob(f"{d}/transfer/Project_*/Data"))
    subdirs.extend(glob.glob(f"{d}/transfer/Project_*/Analysis*"))

    # check if subdir is actually a directory
    for d2 in subdirs:
        if not os.path.isdir(d2):
            subdirs.remove(d2)

    # define txt files to be compressed (legacy)
    compress = glob.glob(f"{d}/*.txt")
    # check if txt file is actually a file
    for f in compress:
        if not os.path.isfile(f):
            compress.remove(f)

    # do nothing if no subdirs
    if len(subdirs) + len(compress) == 0:
        logging.debug(f"{d} has no subdirs to be deleted or removed, skipping")
        continue

    if do_delete:
        logging.info("Tidy " + d)
        for d2 in subdirs:
            logging.info(f"Removing {d2}")
            shutil.rmtree(f"{d2}")
        for f in compress:
            logging.info(f"Compressing {f}")
            subprocess.check_call(["gzip", f])
    else:
        for d2 in subdirs:
            logging.info(f"Would removed: {d2}")
        for f in compress:
            logging.info(f"Would compress {f}")

logging.info("Cleaning done.")
