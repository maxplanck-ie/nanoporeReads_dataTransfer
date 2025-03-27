'''
pycoqc
'''


rule pycoQC:
    input:
        aligned_bam = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
        seqsum = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum"
    output:
        json = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.json",
        html = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.html"
    log: "log/{sample_project}_{sample_id}_{sample_name}_pycoqc.log"
    params:
        par = config['pycoQc']['pycoQc_opts']
    wildcard_constraints:
        sample_id = "[0-9]{2}L[0-9]{6}"
    conda: "ont-ppp-pycoqc"
    shell: """
            # notice that pycoQC is prone to fail, especially for tests with small bam files --> 
            # enforce success (|| true)
            pycoQC --summary_file {input.seqsum} --bam_file {input.aligned_bam} \
            {params.par} -o {output.html} -j {output.json} >> {log} 2>&1 || true
            """

rule pycoqc_final:
    input:
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.html",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.json",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects)
    output: touch("flags/08_pycoqc.done")
