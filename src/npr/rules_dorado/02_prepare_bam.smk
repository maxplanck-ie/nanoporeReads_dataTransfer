import os
import re

idir=config["info_dict"]["base_path"]
baseout=os.path.join(config['info_dict']['flowcell_path'], "bam")

def get_sample_barcode(sample):
    barcode = metadata.loc[metadata['Sample_Name'] == sample, 'index_id'].item()
    return barcode

def get_batches(wildcards):
    # Get checkpoint output (bam_chunks directory)
    checkpoint_output = glob_wildcards(os.path.join(checkpoint_output, "bam_list_{batch}")).batch

    # Extract all batch names (based on chunked files)
    return glob_wildcards(f"{checkpoint_output}/bam_list_*").batch


def collect_batch_bams(wildcards):
    checkpoint_data = checkpoints.split_bam_list_02.get(**wildcards)
    
    checkpoint_output = checkpoint_data.output[0]
    
    batches = glob_wildcards(os.path.join(checkpoint_output, "bam_list_{batch}")).batch
    
    primary_output = expand("bam/{sample_name}_bam_chunks/{batch}.bam", batch=batches, sample_name=wildcards.sample_name)
    #for some reason, snakemake prepends wilcards with whitespaces and results in nonexisting file names
    sanitized_output = [re.sub(r"\s+", "", file) for file in primary_output]
    
    return sanitized_output


rule prepare_bam_list_02:
    input: 
        "flags/01_prepare.done",
        bams_pass = lambda wildcards: os.path.join(idir,"bam_pass") if not demux_done_by_deepseq else os.path.join(idir,"bam_pass", get_sample_barcode(wildcards.sample_name)),
        bams_fail = lambda wildcards: os.path.join(idir,"bam_fail") if not demux_done_by_deepseq else os.path.join(idir,"bam_fail", get_sample_barcode(wildcards.sample_name))
    output:
        bamlist = temp("bam/{sample_name}_bam_list.txt"), #temp
    log: 
        'log/{sample_name}_02_prepare_bam_list.log'
    shell:'''
            # Get bam list from both dirs
            find "{input.bams_pass}" "{input.bams_fail}" -name '*.bam' -type f > "{output.bamlist}" 2>> {log}
        '''


checkpoint split_bam_list_02:
    input:
        bamlist = "bam/{sample_name}_bam_list.txt"
    output:
        chunkdir = temp(directory("bam/{sample_name}_bam_chunks")),
    params:
        batch_size=config['bam_merge']['batch_size']
    log:
        'log/{sample_name}_02_split_bam_list.log'
    shell:'''
            mkdir -p {output.chunkdir}
            # split the list
            split -l "{params.batch_size}" "{input.bamlist}" "{output.chunkdir}/bam_list_b" 2>> {log}
        '''


rule merge_bams_in_batches_02:
    input:
        lambda wildcards: re.sub(r"\s+", "",f"bam/{wildcards.sample_name}_bam_chunks/bam_list_{wildcards.batch}")
    output:
        batch_bam = temp("bam/{sample_name}_bam_chunks/{batch}.bam") #temp
    params:
        opt=config['bam_merge']['opt']
    threads: 40
    conda: "envs/align.yaml"
    log:
        'log/{sample_name}_02_merge_bams_in_batches_{batch}.log'
    benchmark:
        "benchmarks/{sample_name}_00_merge_bams_in_batches_{batch}.tsv"
    shell:
        """
        samtools merge {params.opt} -@ {threads} -b {input} -o {output.batch_bam} 2>> {log}
        """


rule merge_final_bam_02:
    input:
        collect_batch_bams
    output:
        "bam/{sample_id}_{sample_name}.bam"
    params:
        opt=config['bam_merge']['opt']
    threads: 40
    conda: "envs/align.yaml"
    log:
        "log/{sample_id}_{sample_name}_02_merge_final_bam.log"
    shell:
        """
        samtools merge {params.opt} -@ {threads}  -o {output} {input} 2>> {log};
        """


rule prepare_bam_flag_02:
    input: expand("bam/{sample_id}_{sample_name}.bam", zip, sample_id=sample_ids, sample_name=sample_names)
    output: touch("flags/02_prepare_bam.done")
