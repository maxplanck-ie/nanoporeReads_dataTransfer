'''
Extract modification calls after alignment using ONT "modkit"
https://github.com/nanoporetech/modkit
'''

    
rule bam2modbed:
    input:
        flag = "flags/06_align.done",
        bam = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
    output:
        bed = temp("transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed")
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_modbed.log"
    params:
        do_modbed = do_modbed
    threads: 10
    conda:
        "ont-ppp-modkit"
    shell:'''
        echo "do_modbed: {params.do_modbed}" 2>> {log}
        echo modkit pileup -t {threads} {input.bam} {output.bed} 2>> {log}
        modkit pileup -t 100 {input.bam} {output.bed} 2>> {log}
    '''

rule tabix:
    input: "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed"
    output: "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed.gz",
            "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed.gz.tbi"
    log: "log/{sample_project}_{sample_id}_{sample_name}_tabix.log"
    params:
        do_modbed = do_modbed
    conda: "ont-ppp-samtools"
    shell: """
            bgzip {input}
            tabix {output[0]}
        """



rule modbed_final:
    input: expand("transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed.gz.tbi",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/07_modbed.done")
