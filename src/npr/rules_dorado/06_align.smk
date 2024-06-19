'''
Align reads usign dorado aligner (minimap2)
'''
# define source and target pattern
source_bam = sample_dat + ".bam"
source_fq = sample_dat + ".fastq.gz"
source_seqsum = sample_dat + ".seqsum"
target_bam  = sample_ana + ".align.bam"
target_bai  = sample_ana + ".align.bam.bai"
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
        fq_file = source_fq,
        seqsum = source_seqsum
    output:
        file = target_bam,
        bai = target_bai,
        flagstat = target_flagstat
    log:
        logpat
    benchmark:
        bchpat
    params:
        cmd=config['dorado_basecaller']['dorado_cmd'],
        genome =  config['genome'].get(org, None),
        json = target_json,
        html = target_html,
        do_align = config['info_dict']['do_align']
    threads:
        10 # same number of threads will be applied to dorado-align and samtools
    conda:
        "ont-ppp-align"
    shell:
        """
        echo "do_align: {params.do_align}" 2>> {log}
        if [[ "{params.do_align}" == "YES" ]]; then
            {dorado} aligner -t {threads} {params.genome} {input.fq_file} | samtools sort -@ {threads} -m 20G - > {output.file} 2>> {log}
        else
            echo "Alignment step skipped" 2>> {log}
            mv {input.bam_file} {output.file}
        fi

        samtools index -@ {threads} {output.file} 2>> {log}

        # run pycoQC including bam file (for alignment)
        # notice that pycoQC is prone to fail, especially for tests with small bam files --> 
        # enforce success (|| true) and do _not_ require presence of html and json as output of this rule
        pycoQC --summary_file {input.seqsum} --bam_file {output.file} \
        {config[pycoQc][pycoQc_opts]} -o {params.html} -j {params.json} >> {log} 2>&1 || true

        # run samtools flagstat
        samtools flagstat --threads {threads} {output.file} > {output.flagstat}
        """

