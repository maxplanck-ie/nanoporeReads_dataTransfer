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
        idir = config["info_dict"]["base_path"]
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
            echo pod5 merge ${{pod5_files}} -o {output.pod5} 2>> {log}
            pod5 merge ${{pod5_files}} -o {output.pod5} 2>> {log}
        elif [ -e "{params.idir}/fast5_pass" ] || [ -e "{params.idir}/fast5_fail" ]; then
            # there are fast5 produced (legacy)
            fast5_files=`find "{params.idir}/" -name '*.fast5'`
            echo pod5 convert fast5 ${{fast5_files}} --output {output.pod5} 2>> {log}
            pod5 convert fast5 ${{fast5_files}} --output {output.pod5} 2>> {log}
        else
            echo "No raw data found in {params.idir}, hopefully you got BAMs" 2>> {log}
        fi
        '''


