import os
import sys
import glob
import re
import shutil
import pandas as pd
from rich import print

# # define project-specific paths and variable
# org  = config['info_dict']['organism']
# prot = config['info_dict']['protocol']

# create a pandas dataframe of samples to get the sample : project relationship
metadata = dict(config["data"])
del metadata["projects"]
del metadata["samples"]
metadata = pd.DataFrame(metadata).T
sample_ids = metadata['Sample_ID'].tolist()
sample_names = metadata['Sample_Name'].tolist()
sample_projects = metadata['Sample_Project'].tolist()

# Conditional modkit output
expected_modkit = (expand(
        "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.pileup.bed.gz.tbi",
        zip,
        sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects
    ) if (
        config["info_dict"]["modifications"] and # modifications called during basecalling.
        config["info_dict"]["organism_genome"] and # reference genome used during basecalling.
        config["info_dict"]["organism_label"] # Name for the analysis.
    ) else [])

# wildcard mapping from wildcard name -> metadata column
wc_mapping = {
    "project" : "Sample_Project",
    "sample_id" : "Sample_ID",
    "sample_name" : "Sample_Name"
}       

def expand_project_path(path, metadata=metadata, wc_mapping=wc_mapping):
    """
    Given a path containing {project}, {sample_id}, and {sample_name}, 
    find the corresponding columns in the metadata object
    """
    mapping = {k : metadata[v] for k, v in wc_mapping.items()}
    return expand(path, zip, **mapping)
   
# global wildcard constraints: ensure that sample_id adheres to certain constraints: 23L000001
# clarify ambiguities if {sample_id}_{sample_name} = "{23L000001}_{MySample_Part_1}"
wildcard_constraints:
        sample_id="[0-9]{2}L[0-9]{6}",
        sample_name="[a-zA-Z0-9_-]+"

rule finalize:
    input:
        # Explicit output files. Note that these are redundant but kept for clarity (i.e. seqsum implies bam is created, etc.).
        [expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}." + ext, zip, 
               sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects) 
         for ext in ["bam", "bam.bai", "seqsum", "fastq.gz"]],
        "flags/porechop.done",
        expected_modkit,
        # MultiQC (takes fastqc, kraken and pycoqc as inputs.)
        expand(
            "transfer/Project_{sample_project}/QC/multiqc_report.html",
            sample_project = sample_projects
        ),
        # Transfer to periphery.
        "transfer.done"

include: "00_start.smk"
include: "01_prepare.smk"
include: "02_basecall.smk"
include: "03_seqsum.smk"
include: "04_fastq.smk"
include: "05_porechop.smk"
include: "06_modbed.smk"
include: "07_falco.smk"
include: "08_pycoqc.smk"
include: "09_kraken.smk"
include: "10_multiqc.smk"
include: "11_transfer.smk"