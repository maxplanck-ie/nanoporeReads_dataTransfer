'''
Kraken2 contamination test
'''
# define source and target pattern
# source = sample_dat + ".fastq.gz"
# target = sample_qc + "_kraken.report"
# logpat = sample_log + "_kraken.log"
# bchpat = sample_bch + "_kraken.tsv"

rule kraken_final:
    input: expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_kraken.report",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),

    output: touch("flags/08_kraken.done")

rule kraken:
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
    conda:
        "ont-ppp-kraken"
    benchmark:
        "benchmarks/{sample_project}_{sample_id}_{sample_name}_kraken.tsv"
    shell:'''
        kraken2 -db {params.db} --threads {threads} --use-names {input.fastq}  --output {params.output} --report {output.report} 2> {log}
    '''
