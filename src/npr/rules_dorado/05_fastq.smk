'''
Convert BAM files to fastq. Mostly for further use with porechop
Notice that the BAM files may contain modification information
While this can be brute forced into the FASTQ header (-T "*") 
the header will be useless after read chopping, splitting etc
'''

# # define source and target pattern
# source = sample_dat + ".bam" # 
# target = sample_dat + ".fastq.gz" # 
# logpat = sample_log + "_fastq.log"
# bchpat = sample_bch + "_fastq.tsv"

    
rule fastq:
    input:
        flag="flags/04_seqsum.done",
        bam_file = "bam/{sample_id}_{sample_name}.bam"
    output:
        file = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz"
    wildcard_constraints:
        # exclude all sample files that end on ".align.bam" (already aligned) 
        sample_name = r'(?!.*\.align\.bam$).*',
    log:
        "log/{sample_project}_{sample_id}_{sample_name}.fastq.log"
    threads:
        16
    conda:
        "ont-ppp-samtools"
    shell:
        """
        # This assumes that all reads in BAM file are designated READ_OTHER 
        samtools fastq -@ {threads} {input.bam_file} -0 {output.file} 2>> {log}
        """

rule fastq_final:
    input: expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz",zip, sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/05_fastq.done")
