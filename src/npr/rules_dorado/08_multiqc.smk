'''
Conduct various QC-tests
'''

source = ["flags/08_fastqc.done",  "flags/08_kraken.done", "flags/08_pycoqc.done", "flags/05_porechop.done"]
target = os.path.join(project_qc, "multiqc", "multiqc_report.html")
logpat = "log/{project}_multiqc.log"
bchpat = "benchmarks/{project}_multiqc.tsv"

rule multiqc_final:
    input: expand(target, project=config['data']['projects'])
    output: touch("flags/08_multiqc.done")
    shell:'''
        echo "multiqc output: " {input}
    '''

rule multiqc: 
    input: source
    output: 
        target
    params:
        pdir=os.path.join(project_qc, 'Samples'),
        mdir=os.path.join(project_qc, 'multiqc'),
        cfg=config['multiqc']['configfile'],
        sampleDict=os.path.join(project_qc, 'sample_names.tsv'),
        sampleSheet=os.path.join('reports', 'SampleSheet.csv')
    conda:
        "ont-ppp-multiqc"
    log: logpat
    benchmark: bchpat
    shell:'''
        # convert SampleSheet to sampleDict for multiqc convenience
        awk 'BEGIN{{FS=",";OFS="\t"}} !/Data/ {{print $2"_"$3, $2, $3}}' {params.sampleSheet} > {params.sampleDict}

        # do not include 2nd redundant pycoQC report when alignment was run
        multiqc \
            --ignore "*.align_pycoqc*" \
            -c {params.cfg} \
            --sample-names {params.sampleDict} \
            -o {params.mdir} {params.pdir} 2> {log}
    '''