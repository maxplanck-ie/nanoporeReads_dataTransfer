rule fastqc_08: 
    input: 
        "bam/{sample_id}_{sample_name}.bam"
    output: 
        "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_fastqc.html"
    params:
        odir= "transfer/Project_{sample_project}/QC/Samples/",
        memory = config['fastqc'].get('memory', 10000)
    threads:
        10
    conda: "envs/align.yaml"
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_fastqc.log"
    shell:''' 
        fastqc --memory={params.memory} -t {threads} -o {params.odir} --dir {params.odir} {input} 2> {log}
    '''

rule fastqc_final_08: 
    input: expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_fastqc.html",zip, sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects)
    output: touch("flags/08_fastqc.done")
