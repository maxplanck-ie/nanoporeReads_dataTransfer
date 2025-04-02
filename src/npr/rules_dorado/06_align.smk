'''
Align reads usign dorado aligner (minimap2)
'''
# define source and target pattern
source_bam = sample_dat + ".bam"
source_seqsum = sample_dat + ".seqsum"

def is_sorted(sortOrderFile):
    with open(sortOrderFile,'r') as f:
        if "unsorted" in f.readline():
            return False
        else:
            return True


if do_align:   
    rule align:
        input:
            fq_file = "transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.fastq.gz",
        output:
            file = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam"
        log:
            "log/{sample_project}_{sample_id}_{sample_name}_align.log",
        params:
            cmd=config['dorado_basecaller']['dorado_cmd'],
            genome =  config['genome'].get(org, None),
            dir = lambda wildcards: "transfer/Project_" + wildcards.sample_project  + "/" + analysis_name + "/Samples"
        threads:
            10 # same number of threads will be applied to dorado-align and samtools
        conda:
            "ont-ppp-align"
        shell:
            """
            mkdir -p {params.dir};
            {params.cmd} aligner -t {threads} {params.genome} {input.fq_file} | samtools sort -@ {threads} -m 20G - > {output.file} 2>> {log}
        """

else:
    rule grep_sort_order:
        input: "bam/{sample_id}_{sample_name}.bam"
        output: "bam/{sample_id}_{sample_name}.sort_order.txt"
        log: "log/{sample_id}_{sample_name}_sortOrder.log"
        conda: "ont-ppp-samtools"
        shell: """
            samtools view -H {input} | grep SO > {output} 2>{log}
            """

    rule copy_or_sort_aligned_bam:
        input: 
            bam = "bam/{sample_id}_{sample_name}.bam",
            sortorder = "bam/{sample_id}_{sample_name}.sort_order.txt"
        output: 
            file = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam"
        log: "log/{sample_project}_{sample_id}_{sample_name}_copyOrSort.log"
        params:
            dir = lambda wildcards: "transfer/Project_" + wildcards.sample_project  + "/" + analysis_name + "/Samples"
        threads: 10
        conda: "ont-ppp-samtools"
        shell: """
              set -euxo pipefail &&
              gso=$(grep -c -P 'unsorted|unknown' {input.sortorder} || echo "grep failed with exit code $?")
              echo "grep output: $gso" &&
              #touch {output.file} 2>{log} &&
              if (( $gso > 0 ));then
                samtools sort -@ {threads} -m 20G {input.bam} > {output.file} 2>>{log}
              else rsync -av {input.bam} {output.file} 2>>{log}
              fi
              """


rule index_aligned_bam:
    input: "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
    output: "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam.bai"
    log: "log/{sample_project}_{sample_id}_{sample_name}_index.log"
    threads: 10
    conda:
        "ont-ppp-samtools"
    shell: """
        samtools index -@ {threads} {input} 2>> {log}
        """   
    
rule flagstat:
    input:
        aligned_bam = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam",
        bai = "transfer/Project_{sample_project}/" + analysis_name + "/Samples/{sample_id}_{sample_name}.align.bam.bai"
    output:
        flagstat = "transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align.flagstat"
    log: "log/{sample_project}_{sample_id}_{sample_name}_flagstat.log"
    threads: 10
    conda:
        "ont-ppp-samtools"
    shell: """
          samtools flagstat --threads {threads} {input.aligned_bam} > {output.flagstat} 2>{log}
        """

rule align_final:
    input:
         expand("transfer/Project_{sample_project}/" + analysis_name+ "/Samples/{sample_id}_{sample_name}.align.bam.bai",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
         expand("transfer/Project_{sample_project}/QC/Samples/{sample_id}_{sample_name}.align.flagstat",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    output: touch("flags/06_align.done") 
