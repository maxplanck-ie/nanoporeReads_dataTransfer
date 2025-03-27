'''
Align reads usign dorado aligner (minimap2)
'''
# define source and target pattern
source_bam = sample_dat + ".bam"
source_seqsum = sample_dat + ".seqsum"

if do_align:   
    rule align:
        input:
            fq_file = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz",
        output:
            file = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam"
        log:
            "log/{sample_project}_{sample_id}_{sample_name}_align.log",
        params:
            cmd=config['dorado_basecaller']['dorado_cmd'],
            genome =  config['genome'].get(org, None)
        threads:
            10 # same number of threads will be applied to dorado-align and samtools
        conda:
            "ont-ppp-align"
        shell:
            """
            {params.cmd} aligner -t {threads} {params.genome} {input.fq_file} | samtools sort -@ {threads} -m 20G - > {output.file} 2>> {log}
        """

else:
    rule copy_aligned_bam:
        input: "bam/{sample_id}_{sample_name}.bam"
        output: "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam"
        log: "log/{sample_project}_{sample_id}_{sample_name}_rsync.log"
        shell: """
              rsync -av {input} {output} 2>{log}
              """


rule index_aligned_bam:
    input: "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
    output: "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam.bai"
    log: "log/{sample_project}_{sample_id}_{sample_name}_index.log"
    threads: 10
    conda:
        "ont-ppp-samtools"
    shell: """
        samtools index -@ {threads} {input} 2>> {log}
        """   
    
rule flagstat:
    input:
        aligned_bam = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
        bai = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam.bai"
    output:
        flagstat = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align.flagstat"
    log: "log/{sample_project}_{sample_id}_{sample_name}_flagstat.log"
    threads: 10
    conda:
        "ont-ppp-samtools"
    shell: """
          samtools flagstat --threads {threads} {input.aligned_bam} > {output.flagstat} 2>{log}
        """

rule align_final:
    input:
         expand("transfer/Project_{sample_project}/" + analysis_name+ "/Samples/{sample_id}_{sample_name}.align.bam.bai",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
         expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align.flagstat",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/06_align.done") 
