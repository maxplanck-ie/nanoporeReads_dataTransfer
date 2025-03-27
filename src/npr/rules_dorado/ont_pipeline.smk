import os
import sys
import glob
import re
import shutil
import pandas as pd
from rich import print 
from npr.communication import  send_email

dorado=config['dorado_basecaller']['dorado_cmd']

# define flowcell-wide paths
rule_path = config['paths']['rulesPath']
output_bam = os.path.join("bam", "basecalls.bam")
demux_dir  = "bam_demux" 
transfer_dir = "transfer"

# define project-specific paths and variable
org  = config['info_dict']['organism']
prot = config['info_dict']['protocol']
genome = config['genome'].get(org, None)
analysis_name = "_".join(["Analysis", org, prot])

# common project-level path patterns
project_dir = os.path.join(transfer_dir, "Project_{project}")
data_dir = os.path.join(project_dir, "Data")
project_qc = os.path.join(project_dir, "QC")
analysis_dir = os.path.join(project_dir, analysis_name)

# common sample level patterns
# a sample has a (user-defined) name, a (parkour/samplesheet-defined) ID, and an associated project
# here we define workflow-wide patterns for the paths to hold data, analysis, QC, logs and benchmark files
sample_prefix = "{sample_id}_{sample_name}" 
sample_dat = os.path.join(data_dir, sample_prefix)
sample_ana = os.path.join(analysis_dir, sample_prefix)
sample_qc  = os.path.join(project_qc, "Samples", sample_prefix)
sample_log = os.path.join("log", "{project}_" + sample_prefix)
sample_bch = os.path.join("benchmarks", "{project}_" + sample_prefix)

# create a pandas dataframe of samples to get the sample : project relationship
metadata = dict(config["data"])
del metadata["projects"]
del metadata["samples"]
metadata = pd.DataFrame(metadata).T

sample_ids = metadata['Sample_ID'].tolist()
sample_names = metadata['Sample_Name'].tolist()
sample_projects = metadata['Sample_Project'].tolist()
Project_id = sample_projects[0] ## to be used in the fastqc
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


# make alignment conditional on genome being defined and available
do_align = config['info_dict']['do_align']
do_sort = config['info_dict']['do_sort']

#make demultiplexing conditional on barcoding set to true and demultixplexed folders existing on dont_touch_this
barcoding = config['info_dict']['barcoding']
demux_done_by_deepseq = False
if barcoding:
    #get list of all expected demuxed bam folders on dont_touch_this
    check_dirs = [os.path.join(config["info_dict"]["base_path"], "bam_pass" ,x) for x in metadata['index_id'].tolist()]
    if all([os.path.exists(x) for x in check_dirs]):
        demux_done_by_deepseq = True
    else:
        sys.stderr.write("Barcoding set to true but no matching directories found on dont_touch_this.\n")
        exit(1)

#if protocol is cdna, don't call modifications
do_modbed = config['info_dict']['do_modbed']
protocol = config['info_dict']['protocol']
if protocol == "cdna":
    do_modbed = False
    config['info_dict']['do_modbed'] = False
if do_modbed:
    do_modbed_output = expand("transfer/Project_{sample_project}/" + analysis_name+ "/Samples/{sample_id}_{sample_name}.align.bed.gz.tbi",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects).append("flags/07_modbed.done")
else:
    do_modbed_output=[]
   
# global wildcard constraints: ensure that sample_id adheres to certain constraints: 23L000001
# clarify ambiguities if {sample_id}_{sample_name} = "{23L000001}_{MySample_Part_1}"
wildcard_constraints:
        sample_id="[0-9]{2}L[0-9]{6}"


rule finalize:
    input:
        "flags/00_start.done",
        "flags/00_prepare.done",

        "flags/00_prepare_bam.done", 
        expand("bam/{sample_id}_{sample_name}.bam",zip, sample_id=sample_ids, sample_name=sample_names),
        
        "flags/04_seqsum.done", 
        expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        
        "flags/05_fastq.done", 
        expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        
        "flags/05_porechop.done", 
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_porechop.info",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        
        "flags/06_align.done",
        expand("transfer/Project_{sample_project}/" + analysis_name+ "/Samples/{sample_id}_{sample_name}.align.bam",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        
        do_modbed_output,

        "flags/08_fastqc.done",
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_fastqc.html",zip, sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects),
    
        "flags/08_pycoqc.done",
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.html",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.json",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        
        "flags/08_kraken.done",
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_kraken.report",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),

        "flags/08_multiqc.done",
        "transfer/Project_" + Project_id + "/QC/multiqc/multiqc_report.html",

        "flags/09_transfer.done",

        ###"flags/01_basecall.done", 
        ###"flags/02_demux.done",
        ###"flags/03_rename.done",
    

        
    output:    
         touch("flags/XX_snakemake.done")
    log:
        log="log/ont_pipeline.log",

include: "00_start.smk"
include: "00_prepare.smk"
include: "00_prepare_bam.smk"
include: "04_seqsum.smk"
include: "05_fastq.smk"
include: "05_porechop.smk"
include: "06_align.smk"
include: "07_modbed.smk"
include: "08_fastqc.smk"
include: "08_pycoqc.smk"
include: "08_kraken.smk"
include: "08_multiqc.smk"
include: "09_transfer.smk"

#include: "01_basecall.smk"
#include: "02_demux.smk"
#include: "03_rename.smk"



