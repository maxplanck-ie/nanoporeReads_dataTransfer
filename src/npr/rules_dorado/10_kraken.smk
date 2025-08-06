'''
Kraken2 contamination test
'''

rule kraken_final_10:
    input: expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_kraken.report",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),

    output: touch("flags/10_kraken.done")

rule kraken_10:
    input: 
        fastq="transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz"
    output:
        report="transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_kraken.report"
    params:
        db=config['kraken']['db'],
        output="-"  # supress output
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_kraken.log"
    threads:
        10
    conda: "envs/align.yaml"
    benchmark:
        "benchmarks/{sample_project}_{sample_id}_{sample_name}_kraken.tsv"
    shell:'''
        kraken2 -db {params.db} --threads {threads} --use-names {input.fastq}  --output {params.output} --report {output.report} 2> {log}
    '''
