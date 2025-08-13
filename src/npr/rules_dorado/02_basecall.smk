from rich import print
import subprocess as sp
import sys

rule basecall_02:
    input:
        pod5 = 'pod5/output.pod5'
    output:
        bam = 'bam/output.bam'
    threads: 16
    params:
        cmd = config['dorado_basecaller']['dorado_cmd'],
        model_dir = config['dorado_basecaller']['model_directory'],
        model = config['info_dict']['model'],
        reference = f"--reference {config['info_dict']['organism_genome']}" if config['info_dict']['organism_genome'] else "",
        device = config['dorado_basecaller']['device'],
        demux = f"--kit-name {config['info_dict']['barcode_kit']}" if config['info_dict']['barcoding'] else "",
    shell:'''
    {params.cmd} basecaller \
        --models-directory {params.model_dir} \
        {params.reference} {params.demux} \
        -x {params.device} \
        {params.model} {input.pod5} > {output.bam}
    '''

rule basecall_demux_and_sort_02:
    input:
        bam = 'bam/output.bam'
    output:
        bam = expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.bam",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
        bai = expand("transfer/Project_{sample_project}/Data/{sample_id}_{sample_name}.bam.bai",zip, sample_id=sample_ids,sample_name=sample_names, sample_project=sample_projects),
    params:
        cmd = config['dorado_basecaller']['dorado_cmd'],
        demux = f"config[info_dict']['barcode_kit']" if config['info_dict']['barcoding'] else "",
        demuxdir = temp('demux')
    threads: 16
    conda: "envs/align.yaml"
    run:
        # If barcoding is appropriate here, we need to run dorado demux
        # Else, we can just rename the bam files.
        # Note that dorado output is sorted by default, though this additional pass makes it independent of dorado versions.

        if config["info_dict"]["barcode_kit"] == "no_bc":
            assert len(sample_ids) == 1
            assert len(sample_names) == 1
            obam = f"transfer/Project_{sample_projects[0]}/Data/{sample_ids[0]}_{sample_names[0]}.bam"
            shell(f"samtools sort -o {obam} -@ {threads} -m 20G {input.bam}")
            shell(f"samtools index -@ {threads} {obam}")
        else:
            from pathlib import Path
            assert config['info_dict']['barcoding']
            
            shell(f"{params.cmd} demux --no-classify -t {threads} -o {params.demuxdir} {input.bam}")
            # Since this is done, the demux'ed samples will be present in the demuxmdir.
            # We iterate over our metadata table, define output per relevant barcode, and again sort and index
            bamfiles = list(Path(params.demuxdir).glob("*.bam"))
            for i,r in metadata.iterrows():
                bcname = r['index_id']
                bamf = [b for b in bamfiles if bcname in b.name]
                assert len(bamf) == 1, f"Expected one bam file for barcode {bcname}, found {bamf}"
                obam = f"transfer/Project_{r['Sample_Project']}/Data/{r['Sample_ID']}_{r['Sample_Name']}.bam"
                shell(f"samtools sort -o {obam} -@ {threads} -m 20G {bamf[0]}")
                shell(f"samtools index -@ {threads} {obam}")
