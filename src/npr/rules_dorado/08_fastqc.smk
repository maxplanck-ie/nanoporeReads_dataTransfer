'''
Run fastqc
Notice that this runs for all samples in a project (*bam and/or *align.bam)
not sample-by-sample
'''

# define source and target patterns
source = sample_dat + ".bam"
target = sample_qc + "_fastqc.html"
logpat = sample_log + "_fastqc.log"
bchpat = sample_bch + "_fastqc.tsv"
qc_dir = os.path.join(project_dir, "QC", "Samples")

rule fastqc_final: 
    input: expand_project_path(target)
    output: touch("flags/08_fastqc.done")

rule fastqc: 
    input: 
        source
    output: 
        target
    params:
        odir = lambda wildcards: qc_dir.format(project=wildcards.project)
    threads: 4
    log:
        logpat
    benchmark:
        bchpat
    shell:''' 
        fastqc --memory=4096 -t {threads} -o {params.odir} --dir {params.odir} --quiet {input} 2> {log}
    '''