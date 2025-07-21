'''
Prepare BAM file
As the new ONT machine can provide basecalls using dorado + genome
This part checks if we have the BAM files and merge
'''

import os
import re

idir=config["info_dict"]["base_path"]
baseout=os.path.join(config['info_dict']['flowcell_path'], "bam")

def get_sample_barcode(sample):
    barcode = metadata.loc[metadata['Sample_Name'] == sample, 'index_id'].item()
    #print (barcode)
    return barcode

rule prepare_bam_list:
    input: 
        "flags/00_prepare.done",
        bams_pass = lambda wildcards: os.path.join(idir,"bam_pass") if not demux_done_by_deepseq else os.path.join(idir,"bam_pass", get_sample_barcode(wildcards.sample_name)),
        bams_fail = lambda wildcards: os.path.join(idir,"bam_fail") if not demux_done_by_deepseq else os.path.join(idir,"bam_fail", get_sample_barcode(wildcards.sample_name))
    output:
        bamlist = temp("bam/{sample_name}_bam_list.txt"), #temp
    log: 
        'log/{sample_name}_00_prepare_bam_list.log'
    shell:'''
            # Get bam list from both dirs
            find "{input.bams_pass}" "{input.bams_fail}" -name '*.bam' -type f > "{output.bamlist}" 2>> {log}
        '''


checkpoint split_bam_list:
    input:
        bamlist = "bam/{sample_name}_bam_list.txt"
    output:
        chunkdir = temp(directory("bam/{sample_name}_bam_chunks")),
    params:
        batch_size=config['bam_merge']['batch_size']
    log:
        'log/{sample_name}_00_split_bam_list.log'
    shell:'''
            mkdir -p {output.chunkdir}
            # split the list
            split -l "{params.batch_size}" "{input.bamlist}" "{output.chunkdir}/bam_list_b" 2>> {log}
        '''


def get_batches(wildcards):
    # Get checkpoint output (bam_chunks directory)
    checkpoint_output = glob_wildcards(os.path.join(checkpoint_output, "bam_list_{batch}")).batch

    # Extract all batch names (based on chunked files)
    return glob_wildcards(f"{checkpoint_output}/bam_list_*").batch


rule merge_bams_in_batches:
    input:
        lambda wildcards: re.sub(r"\s+", "",f"bam/{wildcards.sample_name}_bam_chunks/bam_list_{wildcards.batch}")
    output:
        batch_bam = temp("bam/{sample_name}_bam_chunks/{batch}.bam") #temp
    params:
        opt=config['bam_merge']['opt']
    threads: 10
    conda:
        "ont-ppp-samtools"
    log:
        'log/{sample_name}_00_merge_bams_in_batches_{batch}.log'
    benchmark:
        "benchmarks/{sample_name}_00_merge_bams_in_batches_{batch}.tsv"
    shell:
        """
        #samtools merge {params.opt} -@ {threads} -b {input} -o {output.batch_bam} 2>> {log}
        samtools merge {params.opt} -@ 100 -b {input} -o {output.batch_bam} 2>> {log}
        """


def collect_batch_bams(wildcards):
    checkpoint_data = checkpoints.split_bam_list.get(**wildcards)
    #print(f"Checkpoint output: {checkpoint_data.output}",file=sys.stderr, flush=True)

    checkpoint_output = checkpoint_data.output[0]
    #print(f"Using checkpoint output: {checkpoint_output}",file=sys.stderr, flush=True)

    batches = glob_wildcards(os.path.join(checkpoint_output, "bam_list_{batch}")).batch
    #print(f"Found batches: {batches}",file=sys.stderr, flush=True)

    primary_output = expand("bam/{sample_name}_bam_chunks/{batch}.bam", batch=batches, sample_name=wildcards.sample_name)
    #for some reason, snakemake prepends wilcards with whitespaces and results in nonexisting file names
    sanitized_output = [re.sub(r"\s+", "", file) for file in primary_output]
    #print(f"Sanitized output: {sanitized_output}",file=sys.stderr, flush=True)
    return sanitized_output


rule merge_final_bam:
    input:
        collect_batch_bams
    output:
        "bam/{sample_id}_{sample_name}.bam"
    params:
        opt=config['bam_merge']['opt']
    threads: 10
    conda:
        "ont-ppp-samtools"
    log:
        "log/{sample_id}_{sample_name}_00_merge_final_bam.log"
    shell:
        """
        samtools merge {params.opt} -@ 100  -o {output} {input} 2>> {log};
        """


rule prepare_bam_flag:
    input: expand("bam/{sample_id}_{sample_name}.bam",zip, sample_id=sample_ids, sample_name=sample_names)
    output: touch("flags/00_prepare_bam.done")
