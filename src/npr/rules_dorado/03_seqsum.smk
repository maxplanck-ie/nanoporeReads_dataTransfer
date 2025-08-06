rule seqsum_03:
    input:
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

rule seqsum_final_03:
    input:
        expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        "flags/02_basecall.done" if do_basecall else "flags/02_prepare_bam.done"
    output:
        touch("flags/03_seqsum.done")
