rule demux:
    input:
        "flags/01_basecall.done"
    output:
        touch("flags/02_demux.done")
    log:
        "log/02_demux.log"
    params:
        kit=config['bc_kit'],
        dir=directory(demux_dir),
        inp=output_bam
    benchmark:
        "benchmarks/02_demux.tsv"
    shell:
        """
        if [ {params.kit} == "no_bc" ]; then
            echo "no demultiplexing" 2>> {log}
        else
            echo "demutliplexing with {params.kit}" 2>> {log}
            {dorado} demux --kit-name {params.kit} --output-dir {params.dir} {params.inp} 2>> {log}        
        fi  
        """

### two unused functions ###
def get_demux_files_notused(wildcards):
    '''
    return a list of all filenames generated by demux (without BAM suffix) - unknown before demux rule is done
    This function is not used any longer. If this functionality is to be revived then change: 
    rule demux --> checkpoint demux
    '''
    ck_output = checkpoints.demux.get(**wildcards).output[0]
    FILES, = glob_wildcards(os.path.join(ck_output, "{filename}.bam"))
    return FILES
    
def seqsum_files_notused(wildcards):
    '''
    return a list of all seqsum files (given a list of demux files)
    use get_demux_files to define various derived filenames (expected output = input for final rule)
    '''
    return expand( seqsum_dir + "/{fn}.seqsum", fn=get_demux_files_notused(wildcards))