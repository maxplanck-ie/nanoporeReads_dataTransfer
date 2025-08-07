rule pycoQC_09:
    input:
        aligned_bam = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
        aligned_bai = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam.bai",
        seqsum = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum"
    output:
        json = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.json",
        html = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.html"
    log: "log/{sample_project}_{sample_id}_{sample_name}_pycoqc.log"
    params:
        par = config['pycoQc']['pycoQc_opts']
    wildcard_constraints:
        sample_id = "[0-9]{2}L[0-9]{6}"
    conda: "envs/pycoqc.yaml"
    shell: """
    pycoQC --summary_file {input.seqsum} --bam_file {input.aligned_bam} \
    {params.par} -o {output.html} -j {output.json} 2>&1 | tee -a {log}
    """

rule pycoqc_final_09:
    input:
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.html",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.json",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects)
    output: touch("flags/09_pycoqc.done")
