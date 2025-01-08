# define target pattern
target = sample_dat + ".bam"  # "transfer/Project_{project}/Sample_{sample_id}/{sample_name}.bam"
logpat = sample_log + "_rename.log"
bchpat = sample_bch + "_rename.tsv"

rule rename_final:
    input: expand_project_path(target),  
    output: touch("flags/03_rename.done")
    shell:'''
        echo "copied: {input}"
    '''

rule rename:
    input:
        flag="flags/02_demux.done"
    output:
        bam=target,
    params:
        dir=demux_dir
    log:
        logpat 
    benchmark:
        bchpat
    run:
        # default file if no barcoding
        barcode_file = output_bam 
        if config['info_dict']['barcoding']:
            # if barcoding is in place, get barcode from sample_id (i.e. from sample sheet)
            ss_dict = config['data'].get(wildcards.sample_id)
            if ss_dict is None:
                print(f'[red] SampleSheet has no entry for sample_id: {wildcards.sample_id}[/red]')
                exit(1)

            barcode = ss_dict.get('index_id', None)
            if barcode is None:
                print(f'[red] Barcode could not be found for sample_id: {wildcards.sample_id}[/red]')
                exit(1)

            barcode_file = os.path.join(
                params.dir, 
                config['bc_kit'] + "_" + barcode + ".bam"
            )
            
        if not os.path.exists(barcode_file):
            print(f'[red] barcode_file does not exist: {barcode_file}. For sample: {wildcards.sample_id} {wildcards.sample_name}[/red]')
            exit(1)
        
        if ' ' in output.bam:
            # Replace empty spaces when PI has 2 last names
            output.bam.replace(' ', '')
        
        print(f'cp {barcode_file} {output.bam}')
        shell("cp {barcode_file} {output.bam}")
