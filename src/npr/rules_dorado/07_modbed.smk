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
        bedgz = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed.gz"
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_modbed.log"
    params:
        bed = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed",
        do_modbed = config['info_dict']['do_modbed'],
        protocol = config['info_dict']['protocol']

    benchmark:
        "benchmarks/{sample_project}_{sample_id}_{sample_name}_modbed.tsv"
    threads: 10
    conda:
        "ont-ppp-modkit"
    shell:'''
        echo "do_modbed: {params.do_modbed}" 2>> {log}
        if [[ "{params.do_modbed}" == "do_modbed" && "{params.protocol}" != "cdna" ]]; then
            echo modkit pileup -t {threads} {input.bam} {params.bed} 2>> {log}
            modkit pileup -t 100 {input.bam} {params.bed} 2>> {log}
            bgzip {params.bed}
            tabix {output.bedgz}
        else
            echo "Modbed step skipped" 2>> {log}
            touch {output.bedgz}
            echo "No modifications observed" >> {output.bedgz}
            
        fi
    '''

rule modbed_final:
    input: expand("transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bed.gz",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/07_modbed.done")
