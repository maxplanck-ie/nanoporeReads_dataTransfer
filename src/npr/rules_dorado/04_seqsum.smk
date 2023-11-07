# define source and target pattern
source = sample_dat + ".bam"
target = sample_dat + ".seqsum" 
logpat = sample_log + "_seqsum.log"
bchpat = sample_bch + "_seqsum.tsv"

rule seqsum_final:
    input: expand_project_path(target)
    output: touch("flags/04_seqsum.done")    
    
rule seqsum:
    input:
        flag="flags/03_rename.done",
        bam_file=source
    output:
        seqsum = target
    params:
        cmd=config['dorado_basecaller']['dorado_cmd']
    log:
        logpat
    benchmark:
        bchpat
    shell:
        """
        {params.cmd} summary {input.bam_file}  > {output.seqsum} 2>> {log}
        """

