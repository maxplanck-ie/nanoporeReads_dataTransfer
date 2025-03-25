'''
Run fastqc
Notice that this runs for all samples in a project (*bam and/or *align.bam)
not sample-by-sample
'''

# # define source and target patterns
# source = sample_dat + ".bam"
# target = sample_qc + "_fastqc.html"
# logpat = sample_log + "_fastqc.log"
# bchpat = sample_bch + "_fastqc.tsv"
qc_dir = os.path.join(project_dir, "QC", "Samples") 

rule fastqc: 
    input: 
        "bam/{sample_id}_{sample_name}.bam"
    output: 
        #"transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_fastqc.html"
        "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_fastqc.html"
    params:
        #odir = lambda wildcards: qc_dir.format(project=wildcards.project),
        odir= "transfer/Project_{sample_project}/QC/Samples/",
        memory = config['fastqc'].get('memory', 10000)
    threads:
        10
    conda:
        "ont-ppp-fastqc"
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_fastqc.log"
    shell:''' 
        fastqc --memory={params.memory} -t {threads} -o {params.odir} --dir {params.odir} {input} 2> {log}
    '''

rule fastqc_final: 
    input: expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_fastqc.html",zip, sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects)
    output: touch("flags/08_fastqc.done")
