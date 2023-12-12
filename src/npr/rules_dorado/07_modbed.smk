'''
Extract modification calls after alignment using "modbam2bed"
'''

# define source and target pattern
source = sample_ana + ".align.bam"
target = sample_ana + ".align.bed.gz"
logpat = sample_log + "_modbed.log"
bchpat = sample_bch + "_modbed.tsv"
    
rule modbed_final:
    input: expand_project_path(target)
    output: touch("flags/07_modbed.done")
    
rule bam2modbed:
    input:
        flag = "flags/06_align.done",
        bam = source
    output:
        bed = target
    conda:
        # mamba create -n modbam2bed -c bioconda -c conda-forge -c epi2melabs modbam2bed
        "envs/modbam2bed.yaml"
        #"/localenv/pipegrp/anaconda/miniconda3/envs/modbam2bed/"
    log:
        logpat
    benchmark:
        bchpat
    params:
        genome = config['genome'].get(org, None)
    threads: 4
    shell:'''
        # future: consider filtering zero modification or "nan" with "| awk '$11 != "nan" && $11>0.0'"
        #   remove or reduce stderr from processing
        ( modbam2bed -t {threads} --combine {params.genome} {input.bam} | pigz -p {threads} > {output.bed} ) 2>> {log}
    '''

