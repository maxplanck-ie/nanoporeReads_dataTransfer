'''
pycoqc
'''

# # define source and target pattern
# source = sample_dat + ".seqsum"
# target_html = sample_qc + "_pycoqc.html"
# target_json = sample_qc + "_pycoqc.json"
# logpat = sample_log + "_pycoqc.log"
# bchpat = sample_bch + "_pycoqc.tsv"

rule pycoqc:
    input: 
        seqsum = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum"
    output:
        html = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.html",
        json = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.json",
    wildcard_constraints:
        # exclude all sample files that end on ".align.bam" (already aligned) 
        sample_name = r'(?!.*\.align_pycoqc\.json$).*',
        sample_id = "[0-9]{2}L[0-9]{6}"
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_pycoqc.log"
    params:
        par = config['pycoQc']['pycoQc_opts']
    conda:
        "ont-ppp-pycoqc"
    shell:'''
        touch {output.html}  # since pycoqc may fail
        touch {output.json}  # since pycoqc may fail
        pycoQC --summary_file {input.seqsum} {params.par} -o {output.html} -j {output.json} >> {log} 2>&1 || true
    '''

rule pycoqc_final:
    input:
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.html",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}_pycoqc.json",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects)
    output: touch("flags/08_pycoqc.done")
