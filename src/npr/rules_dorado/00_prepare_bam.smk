'''
Prepare BAM file
As the new ONT machine can provide basecalls using dorado + genome
This part checks if we have the BAM files and merge
'''

import os

idir=config["info_dict"]["base_path"]
baseout=os.path.join(config['info_dict']['flowcell_path'], "bam")

def get_sample_barcode(sample):
    barcode = metadata['Sample_Name' == sample].'index_id'
    return barcode


rule prepare_bam:
    input: 
        "flags/00_prepare.done",
        bams_pass = lambda wildcards: os.path.join(idir,"bam_pass") if not demux_done_by_deepseq else os.path.join(idir,"bam_pass", get_sample_barcode(wildcards.sample)),
        bams_fail = lambda wildcards: os.path.join(idir,"bam_fail") if not demux_done_by_deepseq else os.path.join(idir,"bam_fail", get_sample_barcode(wildcards.sample))
    output:
        flag=touch("flags/00_prepare_bam.done"),
        bam = os.path.join(baseout,"{sample}_basecalls.bam")
    params:
        baseout=os.path.join(config['info_dict']['flowcell_path'], "bam"),
        batch_size=config['bam_merge']['batch_size'],
        opt=config['bam_merge']['opt']
    threads:
        10
    conda:
        "ont-ppp-samtools"
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

        if [ -e "{input.bams_pass}" ] || [ -e "{input.bams_fail}" ]; then
            # BAM files are present in bam_pass and/or bam_fail

            # Get bam list from both dirs
            echo find "{input.bams_pass}" "{input.bams_fail}" -name '*.bam' -type f > "{params.baseout}/bam_list.txt" 2>> {log}
            find "{input.bam_pass}" "{input.bam_fail}" -name '*.bam' -type f > "{params.baseout}/bam_list.txt" 2>> {log}

            # split the list
            echo split -l {params.batch_size} "{params.baseout}/bam_list.txt" "{params.baseout}/bam_list_b" 2>> {log}
            split -l {params.batch_size} "{params.baseout}/bam_list.txt" "{params.baseout}/bam_list_b" 2>> {log}
            
            #and merge the bams
            if [[ $(ls "{params.baseout}"/bam_list_b* | wc -l) -eq "1" ]]; then
                
                echo samtools merge {params.opt} -@ {threads} -b "{params.baseout}/bam_list_baa" -o "{output.bam}" 2>> {log}
                samtools merge {params.opt} -@ {threads} -b "{params.baseout}/bam_list_baa" -o "{output.bam}" 2>> {log}
            else
                # Merge in batches, one by one to avoid opened files limit
                for BATCH in "{params.baseout}"/bam_list_b*; do
                    echo samtools merge {params.opt} -@ {threads} -b $BATCH -o $BATCH.bam 2>> {log}
                    samtools merge {params.opt} -@ {threads} -b $BATCH -o $BATCH.bam 2>> {log}
                done
                # Final merge
                echo samtools merge {params.opt} -@ {threads} -o "{output.bam}" "{params.baseout}"/bam_list_b*.bam 2>> {log}
                samtools merge {params.opt} -@ {threads} -o "{output.bam}" "{params.baseout}"/bam_list_b*.bam 2>> {log}
            fi
        else
            echo "No BAM files found in {input.bams_pass} or {input_bams_fail}" >> {log}
        fi
        rm -f "{params.baseout}"/bam_list*
        '''
