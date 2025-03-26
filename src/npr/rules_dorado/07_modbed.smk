'''
Extract modification calls after alignment using ONT "modkit"
https://github.com/nanoporetech/modkit
'''

# define source and target pattern
# source = sample_ana + ".align.bam"
# target = sample_ana + ".align.bed.gz"
# logpat = sample_log + "_modbed.log"
# bchpat = sample_bch + "_modbed.tsv"
    
    
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
        if [[ "{params.do_modbed}" == "do_modbed" ]]; then
            echo modkit pileup -t {threads} {input.bam} {output.bed} 2>> {log}
            modkit pileup -t 100 {input.bam} {output.bed} 2>> {log}
        else
            echo "Modbed step skipped" 2>> {log}
            touch {output.bed}
            echo "No modifications observed" >> {output.bed}
            
        fi
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
        if [[ "{params.do_modbed}" == "do_modbed" ]]; then
            bgzip {input}
            tabix {output[0]}
        else
            echo "Modbed step skipped" 2>> {log}
            touch {output[1]}
            echo "No modifications observed" >> {output[0]}

        fi
        """



rule modbed_final:
    input: expand("transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed.gz.tbi",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/07_modbed.done")
