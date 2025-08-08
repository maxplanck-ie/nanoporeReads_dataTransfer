rule seqsum_03:
    input:
        bam = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.bam"
    output:
        seqsum = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum"
    params:
        cmd = config['dorado_basecaller']['dorado_cmd']
    log:
        "log/03_{sample_project}_{sample_id}_{sample_name}.seqsum.log"
    shell: """
    {params.cmd} summary {input.bam} > {output.seqsum} 2> >(tee {log} >&2)
    """