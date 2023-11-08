rule start:
    input:
        "summary/ont_pipeline.dag.jpg",
        "summary/disk_space_start.txt" 
    output: 
        touch("flags/00_start.done")


rule dag:
    output:
        jpg="summary/ont_pipeline.dag.jpg",
    log:
        "log/snakemake.command"
    params:
        sfile=config['snakemake']['snakefile'],
        cfile=config['info_dict']['configFile']  
    shell:'''
    
        echo snakemake --snakefile {params.sfile} --configfile {params.cfile} -c1 --use-conda >> {log}
        snakemake --snakefile {params.sfile} --configfile {params.cfile}  -n --dag --rulegraph | dot -Tjpg -o {output.jpg}
    '''

rule check_disk_space:
    output:
        "summary/disk_space_start.txt"
    params:
        # required disk space in kilobytes
        required_space=100000000
    shell:
        """
        # determine available disk space in current work directory
        avail_space=`df -k . | tail -1 | awk '{{print $4}}'`
        
        # exit if too small
        if [ "$avail_space" -lt "{params.required_space}" ]; then
            echo "Insufficient disk space."
            exit 1
        fi

        # write df output to {log} file
        date > {output}
        df -k . | tail -1 >> {output}
        """
