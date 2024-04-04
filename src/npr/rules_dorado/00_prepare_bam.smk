'''
Prepare BAM file
As the new ONT machine can provide basecalls using dorado + genome
This part checks if we have the BAM files and merge
'''

rule prepare_bam:
    input: 
        "flags/00_prepare.done"
    output:
        flag=touch("flags/00_prepare_bam.done"),
    params:
        idir=config["info_dict"]["base_path"],
        baseout=os.path.join(config['info_dict']['flowcell_path'], "bam"),
        batch_size=config['bam_merge']['batch_size'],
        opt=config['bam_merge']['opt']
    threads:    
        10
    log: 
        'log/00_prepare_bam.log'
    benchmark:
        "benchmarks/00_prepare_bam.tsv"
    shell:'''
        if [ -e "{params.baseout}" ]; then
            echo "{params.baseout} exists, any BAM file will be overwritten" 2>> {log}
        else
            echo mkdir "{params.baseout}" 2>> {log}
            mkdir "{params.baseout}" 2>> {log}
        fi

        if [ -e "{params.idir}/bam_pass" ]; then
            # there are BAM produced
        
            # get list of BAMs to merge
            echo find "{params.idir}/bam_pass" -name '*.bam' "{params.baseout}/bam_list.txt" 2>> {log}
            find "{params.idir}/bam_pass" -name '*.bam' > "{params.baseout}/bam_list.txt"

            # merge BAMs in batches
            echo split -l {params.batch_size} "{params.baseout}/bam_list.txt" "{params.baseout}/bam_list_b" 2>> {log}
            split -l {params.batch_size} "{params.baseout}/bam_list.txt" "{params.baseout}/bam_list_b" 2>> {log}
            if [[ $(ls "{params.baseout}"/bam_list_b* | wc -l) -eq "1" ]]; then
                # only one batch is created, creating final ouput
                echo samtools merge {params.opt} -@ {threads} -b "{params.baseout}/bam_list_baa" -o "{params.baseout}/basecalls.bam" 2>> {log}
                samtools merge {params.opt} -@ {threads} -b "{params.baseout}/bam_list_baa" -o "{params.baseout}/basecalls.bam" 2>> {log}
            
            else
                # merge in batches, one by one to avoid opened files limit
                for BATCH in "{params.baseout}"/bam_list_b*; do
                    echo samtools merge {params.opt} -@ {threads} -b $BATCH -o $BATCH.bam 2>> {log}
                    samtools merge {params.opt} -@ {threads} -b $BATCH -o $BATCH.bam 2>> {log}
                done
                # final merge
                echo samtools merge {params.opt} -@ {threads} -o "{params.baseout}/basecalls.bam" "{params.baseout}"/bam_list_b*.bam 2>> {log}
                samtools merge {params.opt} -@ {threads} -o "{params.baseout}/basecalls.bam" "{params.baseout}"/bam_list_b*.bam 2>> {log}
            fi

            # clean up
            echo rm "{params.baseout}"/bam_list* 2>> {log}
            rm "{params.baseout}"/bam_list* 2>> {log}
        else
            echo "No BAM data found in {params.idir}" 2>> {log}
        fi
        '''
