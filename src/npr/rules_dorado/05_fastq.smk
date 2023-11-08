'''
Convert BAM files to fastq. Mostly for further use with porechop
Notice that the BAM files may contain modification information
While this can be brute forced into the FASTQ header (-T "*") 
the header will be useless after read chopping, splitting etc
'''

# define source and target pattern
source = sample_dat + ".bam" # "transfer/Project_{project}/Sample_{sample_id}/{sample_name}.bam"
target = sample_dat + ".fastq.gz" # "transfer/Project_{project}/Sample_{sample_id}/{sample_name}.fastq.gz"
logpat = sample_log + "_fastq.log"
bchpat = sample_bch + "_fastq.tsv"

rule fastq_final:
    input: expand_project_path(target)
    output: touch("flags/05_fastq.done")
    
rule fastq:
    input:
        flag="flags/03_rename.done",
        bam_file = source
    output:
        file = target
    wildcard_constraints:
        # exclude all sample files that end on ".align.bam" (already aligned) 
        sample_name = r'(?!.*\.align\.bam$).*',
    log:
        logpat
    benchmark:
        bchpat
    threads:
        16
    shell:
        """
        # This assumes that all reads in BAM file are designated READ_OTHER 
        # -T "*" writes flags (e.g modifications) into fastq header
        samtools fastq -@ {threads} -T "*" {input.bam_file} -0 {output.file} 2>> {log}

        #the following is more general and has better compression, but it is much slower: 
        #samtools fastq -@ {threads} {input.bam_file} | pigz -p {threads} > {output.file} 2>> {log}
        """

