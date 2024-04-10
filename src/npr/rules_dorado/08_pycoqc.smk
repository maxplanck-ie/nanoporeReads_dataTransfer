'''
pycoqc
'''

# define source and target pattern
source = sample_dat + ".seqsum"
target_html = sample_qc + "_pycoqc.html"
target_json = sample_qc + "_pycoqc.json"
logpat = sample_log + "_pycoqc.log"
bchpat = sample_bch + "_pycoqc.tsv"

rule pycoqc_final:
    input: expand_project_path(target_html)
    output: touch("flags/08_pycoqc.done")

rule pycoqc:
    input: 
        seqsum = source
    output:
        html = target_html,
        json = target_json
    wildcard_constraints:
        # exclude all sample files that end on ".align.bam" (already aligned) 
        sample_name = r'(?!.*\.align_pycoqc\.json$).*',
        sample_id = "[0-9]{2}L[0-9]{6}"
    log: logpat
    params:
        par = config['pycoQc']['pycoQc_opts']
    benchmark:
        bchpat
    conda:
        "ont-ppp-align"
    shell:'''
        #touch {output.html}  # since pycoqc may fail
        pycoQC --summary_file {input.seqsum} {params.par} -o {output.html} -j {output.json} >> {log} 2>&1
    '''