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
        seqsum=source
    output:
        html=target_html,
    wildcard_constraints:
        # exclude all sample files that end on ".align.bam" (already aligned) 
        sample_name = r'(?!.*\.align_pycoqc\.json$).*',
        sample_id="[0-9]{2}L[0-9]{6}"
    log:
        file=logpat,
        json=target_json
    benchmark:
        bchpat
    shell:'''
        touch {output.html}  # since pycoqc may fail
        pycoQC --summary_file {input.seqsum} \
        {config[pycoQc][pycoQc_opts]} -o {output.html} -j {log.json} >> {log.file} 2>&1 || true
    '''