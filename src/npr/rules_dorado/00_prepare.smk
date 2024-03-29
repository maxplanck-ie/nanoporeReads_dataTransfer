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