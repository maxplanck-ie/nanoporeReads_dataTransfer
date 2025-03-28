rule seqsum:
    input:
        flag="flags/00_prepare_bam.done",
        bam_file= "bam/{sample_id}_{sample_name}.bam"
    output:
        seqsum = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum"
    params:
        cmd=config['dorado_basecaller']['dorado_cmd']
    log:
        "log/{sample_project}_{sample_id}_{sample_name}.seqsum.log"
    shell:
        """
        {params.cmd} summary {input.bam_file}  > {output.seqsum} 2>> {log}
        
        """ 

rule seqsum_final:
    input:
        expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects)
    output:
        touch("flags/04_seqsum.done")
