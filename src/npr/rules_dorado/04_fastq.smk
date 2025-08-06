rule fastq_04:
    input:
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
    conda: "envs/align.yaml"
    shell:
        """
        # This assumes that all reads in BAM file are designated READ_OTHER 
        samtools fastq -@ {threads} {input.bam_file} -0 {output.file} 2>> {log}
        """

rule fastq_final_04:
    input: expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz",zip, sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/04_fastq.done")
