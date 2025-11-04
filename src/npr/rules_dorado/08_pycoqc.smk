rule pycoQC_08:
    input:
        bam = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.bam",
        seqsum = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum"
    output:
        json = "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}.pycoqc.json",
        html = "transfer/Project_{sample_project}/QC/{sample_id}_{sample_name}.pycoqc.html"
    log: "log/{sample_project}_{sample_id}_{sample_name}_pycoqc.log"
    wildcard_constraints:
        sample_id = "[0-9]{2}L[0-9]{6}"
    conda: "envs/pycoqc.yaml"
    run:
        _bamstr = ""        
        if (config['info_dict'].get('organism_genome') and "Nanopore 16S Barcoding Kit" not in config['info_dict']['parkour_protocol']):
            _bamstr = f"--bam_file {input.bam}"

        shell(f"pycoQC --summary_file {input.seqsum} {_bamstr} \
            --min_pass_qual 0 -o {output.html} -j {output.json}"
        )