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
        bam = source,
        do_modbed = config['info_dict']['do_modbed']
    output:
        bed = target
    log:
        logpat
    benchmark:
        bchpat
    threads: 4
    shell:'''
        # future: consider filtering zero modification or "nan" with "| awk '$11 != "nan" && $11>0.0'"
        #   remove or reduce stderr from processing
        if [ "{input.modbed}"  == "true" ]; then
            echo modkit pileup -t {threads} {input.bam} {output.bed} 2>> {log}
            modkit pileup -t {threads} {input.bam} {output.bed} 2>> {log}
        else
            echo "modbed step skip" 2>> {log}
        fi
    '''

