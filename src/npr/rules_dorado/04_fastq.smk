rule fastq_04:
    input:
        bam = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.bam"
    output:
        file = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz"
    log:
        "log/04_{sample_project}_{sample_id}_{sample_name}.fastq.log"
    threads: 16
    conda: "envs/align.yaml"
    shell:"""
    samtools fastq -@ {threads} {input.bam} -0 {output.file} 2>&1 | tee {log}
    """