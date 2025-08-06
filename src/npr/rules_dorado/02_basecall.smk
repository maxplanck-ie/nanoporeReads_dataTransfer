from rich import print
import subprocess as sp

def gpu_available():
    # Check if GPU is available
    prog = "nvidia-smi"
    query = "--query-gpu=count,memory.free,utilization.gpu,utilization.memory"
    monitor = [prog,  query, "--format=csv,noheader,nounits"]
    #monitor = ["python", "-c", "import sys; sys.exit(6)"]     # fake unavailable GPU

    if shutil.which(prog) is None:
        # if nvidia-smi is not available assume that there is no GPU on the system
        # in such cases do not exit. if device cuda* is specified such jobs would fail latter
        print("[yellow]{} is not available. Do not monitor GPU[/yellow]".format(prog)  )
        return True

    res = sp.run(monitor, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    print("GPU status: stdout={} stderr={}\n".format(res.stdout, res.stderr))
    if res.returncode == 0:
        res = res.stdout.strip().split(',') # split csv output
        res = [int(i) for i in res]         # convert to integer
        # consider GPU available if memory.free>10 GB, and utilization.gpu<80% and utilization.memory<80%
        if res[1]>10000 and res[2]<80 and res[3]<80: 
            return True
    else:
        print("[red]{} caused an error while checking GPU status: {} [/red]".format(" ".join(monitor),res.returncode))
        #sys.exit(1)

    print("[yellow]GPU is busy. {} [/yellow]".format(res))
    # consider GPU as not available
    return False

def get_model(prot_type):    
    model = (
        config['info_dict']['model']
        if 'dna' in prot_type 
        else '/localenv/pipegrp/dorado-1.0.2-linux-x64/model/rna004_130bps_hac@v5.0.0'
        if 'rna' in prot_type 
        else None
    )

    return model

def get_mod(prot_type):
    mod = (
        config['dorado_basecaller'].get("dorado_dna_modifications", "")
        if "dna" in prot_type 
        else config['dorado_basecaller'].get("dorado_rna_modifications", "")
        if "rna" in prot_type 
        else ""
    )
    
    return mod

def get_genome():
    org=config['info_dict']['organism']
    ref=config['genome'][org]

    return ref
        

rule basecall_02:
    input:
        flag="flags/01_prepare.done",
    output:
        "bam/{sample_id}_{sample_name}.bam"
    log:
        "log/{sample_id}_{sample_name}_00_merge_final_bam.log"
    params:
        cmd=config['dorado_basecaller']['dorado_cmd'],
        prot_type=config['info_dict']['protocol'],
        model = get_model(config['info_dict']['protocol']),  # FIX: removed extra comma
        mod = get_mod(config['info_dict']['protocol']),
        options=config['dorado_basecaller']['dorado_options'],
        dir='pod5',
        bam=os.path.join(config['info_dict']['flowcell_path'], "bam/{sample_id}_{sample_name}.bam"),
        bam_dir=os.path.join(config['info_dict']['flowcell_path'],"bam"),
        reference=get_genome(),
        do_basecall=config['info_dict']['do_basecall']

    run:
        # now we wait for the GPU only when basecalling is needed
        if config['info_dict']['do_basecall'] == "do_basecall":
            while not gpu_available():
                print("GPU is unavailable - sleep")
                time.sleep(60)

        shell(
        """
        echo "do_basecall: {params.do_basecall}" 2>> {log}
        if [[ "{params.prot_type}" == "rna" ]]; then
            echo "Processing RNA modification" 2>> {log}
            echo {params.cmd} basecaller sup {params.dir} {params.options} {params.mod} --reference {params.reference}  > {params.bam} 2>> {log}
            {params.cmd} basecaller sup {params.dir} {params.options} {params.mod} --reference {params.reference} > {params.bam} 2>> {log}
        else
            echo "Processing DNA modification" 2>> {log}
            echo {params.cmd} basecaller sup {params.dir} {params.options} {params.mod} --reference {params.reference}  > {params.bam} 2>> {log}
            {params.cmd} basecaller sup {params.dir} {params.options} {params.mod} --reference {params.reference} > {params.bam} 2>> {log}
        fi

        """
        )

# rule prepare_bam:
#     output:
#         touch("flags/00_prepare_bam.done")
#     input:
#         # Add any required input files here
#     shell:
#         """
#         # Commands to prepare the BAM file
#         touch {output}
#         """

rule basecall_flag_02:
    input: expand("bam/{sample_id}_{sample_name}.bam", zip, sample_id=sample_ids, sample_name=sample_names)
    output: touch("flags/02_basecall.done")