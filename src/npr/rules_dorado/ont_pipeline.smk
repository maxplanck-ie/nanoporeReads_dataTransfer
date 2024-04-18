import os
import sys
import glob
import re
import shutil
import pandas as pd
from rich import print

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
align_done = []
if genome is not None and os.path.exists(genome):
    align_done = ["flags/06_align.done","flags/07_modbed.done"]
else:
    sys.stderr.write(f"No genome for organism {org}. No alignment will be done\n")

# global wildcard constraints: ensure that sample_id adheres to certain constraints: 23L000001
# clarify ambiguities if {sample_id}_{sample_name} = "{23L000001}_{MySample_Part_1}"
wildcard_constraints:
        sample_id="[0-9]{2}L[0-9]{6}"

rule finalize:
    input:
        "flags/00_start.done",
        "flags/00_prepare.done",
        "flags/00_prepare_bam.done",
        "flags/01_basecall.done",
        "flags/02_demux.done",
        "flags/03_rename.done",
        "flags/04_seqsum.done",
        "flags/05_fastq.done",
        "flags/05_porechop.done",
        align_done,
        "flags/08_fastqc.done",
        "flags/08_pycoqc.done",
        "flags/08_kraken.done",
        "flags/08_multiqc.done",
        "flags/09_transfer.done",
    output:    
        touch("flags/XX_snakemake.done")
    benchmark:
        "benchmarks/ont_pipeline.tsv"
    log:
        log="log/ont_pipeline.log",
 
    params:
        bench_comb = "benchmarks_combined.tsv"
    shell:'''
        # combine all benchmark files
        print_header=true
        for file in benchmarks/*.tsv; do
            filename=$(basename "$file")
            if [ "$print_header" = true ] ; then
                head -n 1 "$file" | sed "s/^/filename\\t/"
                print_header=false
            fi
            tail -n +2 "$file" | awk -v filename="$filename" 'BEGIN {{OFS="\\t"}} {{print filename, $_}}' >> {params.bench_comb}
        done >> {params.bench_comb}
    '''

include: "00_start.smk"
include: "00_prepare.smk"
include: "00_prepare_bam.smk"
include: "01_basecall.smk"
include: "02_demux.smk"
include: "03_rename.smk"
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