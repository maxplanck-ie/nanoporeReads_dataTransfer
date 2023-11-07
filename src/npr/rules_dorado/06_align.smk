'''
Align reads usign dorado aligner (minimap2)
'''
# define source and target pattern
source_bam = sample_dat + ".bam"
source_seqsum = sample_dat + ".seqsum"
target_bam  = sample_ana + ".align.bam"
target_html = sample_qc + ".align_pycoqc.html"
target_json = sample_qc + ".align_pycoqc.json"
target_flagstat = sample_qc + ".align.flagstat"
logpat = sample_log + "_align.log"
bchpat = sample_bch + "_align.tsv"

rule align_final:
    input: expand_project_path(target_bam)
    output: touch("flags/06_align.done")
    
rule align:
    input:
        bam_file = source_bam,
        seqsum = source_seqsum
    output:
        file = target_bam,
#        html = target_html,
#        json = target_json,
        flagstat = target_flagstat
#    wildcard_constraints:
#        # exclude all sample files that end on ".align.bam" (already aligned) 
#        sample_name = r'(?!.*\.align\.bam$).*'
    log:
        logpat
    benchmark:
        bchpat
    params:
        cmd=config['dorado_basecaller']['dorado_cmd'],
        genome =  config['genome'].get(org, None),
        json = target_json,
        html = target_html,
    threads:
        4
    shell:
        """
        {dorado} aligner {params.genome} {input.bam_file} | samtools sort - > {output.file} 2>> {log}
        samtools index {output.file} 2>> {log}

        # run pycoQC including bam file (for alignment)
        # notice that pycoQC is prone to fail, especially for tests with small bam files --> 
        # enforce success (|| true) and do _not_ require presence of html and json as output of this rule
        pycoQC --summary_file {input.seqsum} --bam_file {output.file} \
        {config[pycoQc][pycoQc_opts]} -o {params.html} -j {params.json} >> {log} 2>&1 || true

        # run samtools flagstat
        samtools flagstat --threads {threads} {output.file} > {output.flagstat}
        """
