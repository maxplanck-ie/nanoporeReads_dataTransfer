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
