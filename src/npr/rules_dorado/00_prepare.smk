'''
Prepare pod5 file
Starting point of every workflow should be a single pod5 file: pod5/output.pod5
The latter can be merged from all pod5 files at origin = input directory "idir"
1. from pod5: pod5 merge ... (specify output file explicitly)
2. from fast5: pod5-convert-from-fast5 ... (specify output directory )
'''

rule prepare:
    input: 
        "flags/00_start.done"
    output:
        pod5="pod5/output.pod5", 
        flag=touch("flags/00_prepare.done")
    params:
        idir = config["info_dict"]["base_path"],
        baseout = os.path.join(config['info_dict']['flowcell_path'],"pod5"),
        opt = '--recursive --force-overwrite'
    threads:    
        10
    log: 
        'log/00_prepare.log'
    benchmark:
        "benchmarks/00_prepare.tsv"
    shell:'''
        if [ -e "{params.idir}/pod5_pass" ] || [ -e "{params.idir}/pod5_fail" ]; then
            # there are pod5 produced (default)
            pod5_files=`find "{params.idir}/" -name '*.pod5'`
            echo pod5 merge ${{pod5_files}} {output.pod5} 2>> {log}
            pod5 merge ${{pod5_files}} {output.pod5} 2>> {log}
        elif [ -e "{params.idir}/fast5_pass" ] || [ -e "{params.idir}/fast5_fail" ]; then
            # there are fast5 produced (legacy)
            pod5-convert-from-fast5 {params.opt} {params.idir} {params.baseout}  -p {threads} > {log}
        else
            echo "No raw data found in {params.idir}" 2>> {log}
            exit
        fi
        '''

rule prepare_bam:
    input: 
        "flags/00_start.done"
    output:
        dflag=touch("flags/00_prepare_bam.done"),
        flag=touch("flags/01_basecall.done") if os.path.exists(os.path.join("{params.baseout}", "basecalls.bam")) else []
    params:
        idir = config["info_dict"]["base_path"],
        baseout = os.path.join(config['info_dict']['flowcell_path'],"bam"),
        batch_size = 500,
        opt = '-c --no-PG -@'
    threads:    
        10
    log: 
        'log/00_prepare_bam.log'
    benchmark:
        "benchmarks/00_prepare_bam.tsv"
    shell:'''
        if [ -e "{params.idir}/bam_pass" ]; then
            # there are BAM produced (default)
            if [ -e "{params.baseout}" ]; then
                echo "{params.baseout} exists, any BAM file will be overwritten" 2>> {log}
            else
                echo mkdir "{params.baseout}" 2>> {log}
                mkdir "{params.baseout}" 2>> {log}
            fi

            # get list of BAMs to merge
            echo find "{params.idir}/bam_pass" -name '*.bam' > "{params.baseout}/bam_list.txt" 2>> {log}
            find "{params.idir}/bam_pass" -name '*.bam' > "{params.baseout}/bam_list.txt"
            # merge BAMs in batches
            echo split -l {params.batch_size} "{params.baseout}/bam_list.txt" "{params.baseout}/bam_list_b" 2>> {log}
            split -l {params.batch_size} "{params.baseout}/bam_list.txt" "{params.baseout}/bam_list_b" 2>> {log}
            for BATCH in "{params.baseout}"/bam_list_b*; do
                echo samtools merge {params.opt} {threads} -b $BATCH -o $BATCH.bam 2>> {log}
                samtools merge {params.opt} {threads} -b $BATCH -o $BATCH.bam 2>> {log}
            done
            
            # final merge
            if [[ $(ls "{params.baseout}"/bam_list_b*.bam | wc -l) -eq "1" ]]; then
                # only one file was produced, just renaming it
                echo mv "{params.baseout}"/bam_list_b*.bam "{params.baseout}/basecalls.bam" 2>> {log}
                mv "{params.baseout}"/bam_list_b*.bam "{params.baseout}/basecalls.bam" 2>> {log}
            else
                echo samtools merge {params.opt} {threads} -o "{params.baseout}/basecalls.bam" "{params.baseout}"/bam_list_b*.bam 2>> {log}
                samtools merge {params.opt} {threads} -o "{params.baseout}/basecalls.bam" "{params.baseout}"/bam_list_b*.bam 2>> {log}
            fi

            # clean up
            echo rm "{params.baseout}"/bam_list* 2>> {log}
            rm "{params.baseout}"/bam_list* 2>> {log}
        else
            echo "No BAM data found in {params.idir}" 2>> {log}
        fi
        '''
