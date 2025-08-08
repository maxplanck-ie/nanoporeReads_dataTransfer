rule kraken_09:
    input: 
        fastq="transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz"
    output:
        report="transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}_kraken.report"
    params:
        db=config['kraken']['db'],
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_kraken.log"
    threads: 20
    conda: "envs/align.yaml"
    benchmark:
        "benchmarks/{sample_project}_{sample_id}_{sample_name}_kraken.tsv"
    shell:'''
    kraken2 -db {params.db} --threads {threads} --output - --report {output.report} {input.fastq} 2>&1 {log}
    '''
