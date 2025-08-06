rule prepare_01:
    input: 
        "flags/00_start.done"
    output:
        pod5="pod5/output.pod5", 
        flag=touch("flags/01_prepare.done")
    params:
        idir = config["info_dict"]["base_path"]
    threads:
        10
    conda:
        "envs/align.yaml"
    log: 
        'log/01_prepare.log'
    shell:'''
        if [ -e "{params.idir}/pod5_pass" ] || [ -e "{params.idir}/pod5_fail" ] || [ -e "{params.idir}/pod5" ]; then
            # there are pod5 produced (default)
            pod5_files=`find "{params.idir}/" -name '*.pod5'`
            echo pod5 merge ${{pod5_files}} -o {output.pod5} -t {threads} 2>> {log}
            pod5 merge ${{pod5_files}} -o {output.pod5} -t {threads} 2>> {log}
        elif [ -e "{params.idir}/fast5_pass" ] || [ -e "{params.idir}/fast5_fail" ]; then
            # there are fast5 produced (legacy)
            fast5_files=`find "{params.idir}/" -name '*.fast5'`
            echo pod5 convert fast5 ${{fast5_files}} -o {output.pod5} -t {threads} --strict 2>> {log}
            pod5 convert fast5 ${{fast5_files}} -o {output.pod5} -t {threads} --strict 2>> {log}
        else
            echo "No raw data found in {params.idir}, hopefully you got BAMs" 2>> {log}
        fi
        '''