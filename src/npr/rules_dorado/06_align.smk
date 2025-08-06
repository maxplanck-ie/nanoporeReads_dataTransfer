if do_align:   
    rule align_06:
        input:
            fq_file = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz",
        output:
            file = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam"
        log:
            "log/{sample_project}_{sample_id}_{sample_name}_align.log",
        params:
            cmd=config['dorado_basecaller']['dorado_cmd'],
            genome =  config['genome'].get(org, None),
            dir = lambda wildcards: "transfer/Project_" + wildcards.sample_project  + "/" + analysis_name + "/Samples"
        threads:
            10 # same number of threads will be applied to dorado-align and samtools
        conda: "envs/align.yaml"
        shell:"""
        mkdir -p {params.dir};
        {params.cmd} aligner -t {threads} {params.genome} {input.fq_file} | samtools sort -@ {threads} -m 20G - > {output.file} 2>> {log}
        """

else:
    rule sort_aligned_bam_06:
        input: 
            bam = "bam/{sample_id}_{sample_name}.bam",
        output: 
            file = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam"
        log: "log/{sample_project}_{sample_id}_{sample_name}_copyOrSort.log"
        threads: 10
        conda: "envs/align.yaml"
        shell: """
        samtools sort -o {output.file} -@ {threads} -m 20G {input.bam} 2>>{log}
        """

rule index_aligned_bam_06:
    input: 
        "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
    output: 
        "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam.bai"
    log: 
        "log/{sample_project}_{sample_id}_{sample_name}_index.log"
    threads: 10
    conda: "envs/align.yaml"
    shell: """
    samtools index -@ {threads} {input} 2>> {log}
    """   
    
rule flagstat_06:
    input:
        aligned_bam = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
        bai = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam.bai"
    output:
        flagstat = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align.flagstat"
    log: "log/{sample_project}_{sample_id}_{sample_name}_flagstat.log"
    threads: 10
    conda: "envs/align.yaml"
    shell: """
    samtools flagstat --threads {threads} {input.aligned_bam} > {output.flagstat} 2>{log}
    """

rule align_final_06:
    input:
        expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align.flagstat",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    output: 
        touch("flags/06_align.done") 
