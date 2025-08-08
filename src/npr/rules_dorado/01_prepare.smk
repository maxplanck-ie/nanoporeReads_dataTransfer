rule prepare_01:
    input: 
        "summary/00_disk_space_start.txt"
    output:
        pod5="pod5/output.pod5",
    params:
        idir = config["info_dict"]["base_path"]
    threads: 20
    conda: "envs/align.yaml"
    log: 'log/01_prepare.log'
    shell:'''
    pod5_files=`find "{params.idir}/" -name '*.pod5'`
    pod5 merge ${{pod5_files}} -o {output.pod5} -t {threads} 2>&1 | tee {log}
    '''