rule bam2modbed_06:
    input:
        bam = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.bam"
    output:
        bed = temp("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.pileup.bed"),
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_modbed.log"
    threads: 40
    conda: "envs/align.yaml"
    shell:'''
    modkit pileup -t 35 {input.bam} {output.bed}
    '''

rule modbed_compress_06:
    input:
        bed = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.pileup.bed"
    output:
        bed = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.pileup.bed.gz"
    conda: "envs/align.yaml"
    threads: 20
    shell:'''
    bgzip -@ {threads} {input.bed}
    '''

rule modbed_tabix_06:
    input:
        bed = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.pileup.bed.gz"
    output:
        bed = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.pileup.bed.gz.tbi"
    conda: "envs/align.yaml"
    threads: 20
    shell:'''
    tabix -@ {threads} -p bed {input.bed}
    '''