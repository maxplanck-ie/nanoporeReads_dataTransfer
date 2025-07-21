'''
Conduct various QC-tests
'''
source = ["flags/08_fastqc.done",  "flags/08_kraken.done", "flags/08_pycoqc.done", "flags/05_porechop.done"]

rule multiqc_final:
    input: expand("transfer/Project_" + Project_id + "/QC/multiqc/multiqc_report.html")
    output: touch("flags/08_multiqc.done")
    shell:'''
        echo "multiqc output: " {input}
    '''

rule multiqc: 
    input: source
    output:
        "transfer/Project_" + Project_id + "/QC/multiqc/multiqc_report.html"
    params:
        pdir="transfer/Project_" + Project_id + "/QC/Samples/",
        mdir="transfer/Project_" + Project_id + "/QC/multiqc/",
        cfg=config['multiqc']['configfile'],
        sampleDict="transfer/Project_" + Project_id + "/QC/sample_names.tsv",
        sampleSheet="reports/SampleSheet.csv"
    conda:
        "ont-ppp-multiqc"
    log: "log/"+ Project_id+ "_multiqc.log"
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
