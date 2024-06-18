'''
Extract modification calls after alignment using ONT "modkit"
https://github.com/nanoporetech/modkit
'''

# define source and target pattern
source = sample_ana + ".align.bam"
target = sample_ana + ".align.bed.gz"
logpat = sample_log + "_modbed.log"
bchpat = sample_bch + "_modbed.tsv"
    
rule modbed_final:
    input: expand_project_path(target)
    output: touch("flags/07_modbed.done")
    
rule bam2modbed:
    input:
        flag = "flags/06_align.done",
        bam = source
    output:
        bed = target
    log:
        logpat
    params:
        do_modbed = config['info_dict']['do_modbed']
    benchmark:
        bchpat
    threads: 10
    conda:
        "ont-ppp-modkit"
    shell:'''
        echo "do_modbed: {params.do_modbed}" 2>> {log}
        do_modbed={params.do_modbed}
        if [[ ${do_modbed} == "YES" ]]; then
            echo modkit pileup -t {threads} {input.bam} {output.bed} 2>> {log}
            modkit pileup -t {threads} {input.bam} {output.bed} 2>> {log}
        else
            echo "Modbed step skipped" 2>> {log}
        fi
    '''

