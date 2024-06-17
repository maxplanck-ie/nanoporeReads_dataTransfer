'''
Kraken2 contamination test
'''
# define source and target pattern
source = sample_dat + ".fastq.gz"
target = sample_qc + "_kraken.report"
logpat = sample_log + "_kraken.log"
bchpat = sample_bch + "_kraken.tsv"

rule kraken_final:
    input: expand_project_path(target)
    output: touch("flags/08_kraken.done")

rule kraken:
    input: 
        fastq=source
    output:
        report=target
    params:
        db=config['kraken']['db'],
        output="-"  # supress output
    log:
        logpat
    threads:
        10
    conda:
        "ont-ppp-kraken"
    benchmark:
        bchpat
    shell:'''
        kraken2 -db {params.db} --threads {threads} --use-names {input.fastq}  --output {params.output} --report {output.report} 2> {log}
    '''
