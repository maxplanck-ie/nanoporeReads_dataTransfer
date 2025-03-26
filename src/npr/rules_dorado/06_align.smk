'''
Align reads usign dorado aligner (minimap2)
'''
# define source and target pattern
source_bam = sample_dat + ".bam"
# source_fq = sample_dat + ".fastq.gz"
source_seqsum = sample_dat + ".seqsum"
# target_bam  = sample_ana + ".align.bam"
# target_bai  = sample_ana + ".align.bam.bai"
# target_html = sample_qc + ".align_pycoqc.html"
# target_json = sample_qc + ".align_pycoqc.json"
# target_flagstat = sample_qc + ".align.flagstat"
# logpat = sample_log + "_align.log"
# bchpat = sample_bch + "_align.tsv"

    
rule align:
    input:
        bam_file = "bam/{sample_id}_{sample_name}.bam",
        fq_file = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz",
        seqsum = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.seqsum"
    output:
        file = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
        #file = "transfer/Project_{sample_project}/Samples/{sample_id}_{sample_name}.align.bam",
        #bai = "transfer/Project_{sample_project}/{analysis_name}/Samples/{sample_id}_{sample_name}.align.bam.bai",
        flagstat = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align.flagstat"
    log:
        "log/{sample_project}_{sample_id}_{sample_name}_align.log",
    params:
        cmd=config['dorado_basecaller']['dorado_cmd'],
        genome =  config['genome'].get(org, None),
        json = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.json",
        html = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align_pycoqc.html",
        do_align = config['info_dict']['do_align']
    threads:
        10 # same number of threads will be applied to dorado-align and samtools
    conda:
        "ont-ppp-align"
    shell:
        """
        echo "do_align: {params.do_align}" 2>> {log}
        if [[ "{params.do_align}" == "do_align" ]]; then
            {dorado} aligner -t {threads} {params.genome} {input.fq_file} | samtools sort -@ {threads} -m 20G - > {output.file} 2>> {log}
        else
            echo "Alignment step skipped" 2>> {log}
            #mv {input.bam_file} {output.file}
            #sometimes the file is not sorted, add functionality to sort the input file first before indexing
            samtools sort -@ {threads} -m 20G -o {output.file} {input.bam_file}
        fi

        echo "indexing bam {output.file}"
        samtools index -@ {threads} {output.file} 2>> {log}

        echo "pycoQC for bam {output.file}"
        # run pycoQC including bam file (for alignment)
        # notice that pycoQC is prone to fail, especially for tests with small bam files --> 
        # enforce success (|| true) and do _not_ require presence of html and json as output of this rule
        pycoQC --summary_file {input.seqsum} --bam_file {output.file} \
        {config[pycoQc][pycoQc_opts]} -o {params.html} -j {params.json} >> {log} 2>&1 || true

        # run samtools flagstat
        echo "samtools flagstat for bam {output.file}"
        samtools flagstat --threads {threads} {output.file} > {output.flagstat}
        
        """

rule align_final:
    input: expand("transfer/Project_{sample_project}/" + analysis_name+ "/Samples/{sample_id}_{sample_name}.align.bam",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    #input: expand("transfer/Project_{sample_project}/Samples/{sample_id}_{sample_name}.align.bam",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/06_align.done") 
