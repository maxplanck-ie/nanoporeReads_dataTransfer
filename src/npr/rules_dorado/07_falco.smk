rule falco_07: 
    input: 
        "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz"
    output: 
        qcrep = "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}_falqoqc.html",
        qcdata = "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}_fastqc_data.txt",
        qcsum = "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}_fastqc_summary.txt"
    threads:
        20
    conda: "envs/align.yaml"
    shell:'''
    falco -t {threads} -R {output.qcrep} -D {output.qcdata} -S {output.qcsum} {input}
    '''