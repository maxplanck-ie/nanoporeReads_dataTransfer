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
        
rule basecall:
    input:
        flag="flags/00_prepare_bam.done",
    output:
        bam=output_bam,
        flag=touch("flags/01_basecall.done")
    log:
        "log/01_basecall.log"
    benchmark:
        "benchmarks/01_basecall.tsv"
    params:
        cmd=config['dorado_basecaller']['dorado_cmd'],
        model=config['info_dict']['model'],
        options=config['dorado_basecaller']['dorado_options'],
        # modification do not yet work with RNA
        mod=config['dorado_basecaller']['dorado_modifications'] \
            if not config['info_dict']['model_def'].startswith("rna") else "",
        do_basecall=config['info_dict']['do_basecall'],
        dir='pod5'
    run: 
        while not gpu_available():
            print("GPU is unavailable - sleep")
            time.sleep(60)


        shell(
        """
        echo "do_basecall: {params.do_basecall}" 2>> {log}
        if [[ {params.do_basecall} ]]; then
            echo {params.cmd} basecaller {params.model} {params.dir} {params.options} {params.mod} > {output.bam} 2>> {log}
            {params.cmd} basecaller {params.model} {params.dir} {params.options} {params.mod} > {output.bam} 2>> {log}
        else
            echo "Basecall step skip" 2>> {log}
        fi
        """
        )
