rule multiqc_10: 
    input:
        fqc = expand(
            "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}_falqoqc.html",
            sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects
        ),
        pycoqc = expand(
            "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}.pycoqc.html",
            sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects
        ),
        kraken = expand(
            "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}_kraken.report",
            sample_id=sample_ids, sample_name=sample_names, sample_project=sample_projects
        ),
        sampleSheet = "reports/SampleSheet.csv",
        modkit = expected_modkit
    output:
        sampleDict = "transfer/Project_{sample_project}/QC/sample_names.tsv",
        qcrep = "transfer/Project_{sample_project}/QC/multiqc_report.html",
    params:
        pdir= lambda wildcards: f"transfer/Project_{wildcards.sample_project}",
        odir = lambda wildcards: f"transfer/Project_{wildcards.sample_project}/QC",
    conda: "envs/align.yaml"
    shell:'''
    # convert SampleSheet to sampleDict for multiqc convenience
    awk 'BEGIN{{FS=",";OFS="\t"}} !/Data/ {{print $2"_"$3, $2, $3}}' {input.sampleSheet} > {output.sampleDict}
    multiqc \
        --sample-names {output.sampleDict} \
        -o {params.odir} {params.pdir}
    '''
