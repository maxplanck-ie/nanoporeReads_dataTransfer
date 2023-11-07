from npr.snakehelper import guppy_basecalling, dorado_basecalling
from rich import print
import subprocess as sp
import time

def gpu_available():
    # Check if GPU is available
    prog="nvidia-smi"
    query=" --query-gpu=count,memory.free,utilization.gpu,utilization.memory"
    cmd = prog + query + " --format=csv,noheader,nounits"
    if shutil.which(prog) is None:
        # if nvidia-smi is not available assume that there is no GPU on the system
        # in such cases do not exit. if device cuda* is specified such jobs would fail latter
        print("[yellow]{} is not available. Do not monitor GPU[/yellow]".format(prog)  )
        return True

    res = sp.run(cmd.split(), stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    print("GPU status: {}\n".format(res.stdout))
    if res.returncode == 0:
        res = res.stdout.strip().split(',') # split csv output
        res = [int(i) for i in res]         # convert to integer
        # consider GPU available if memory.free>20 MB, and utilization.gpu<80% and utilization.memory<80%
        if res[1]>20000 and res[2]<80 and res[3]<80: 
            return True
    else:
        print("[red]{} caused an error while checking GPU status: {} [/red]".format(cmd,res.returncode))
        #sys.exit(1)

    print("[yellow]GPU is busy. {} [/yellow]".format(res))
    # consider GPU as not available
    return False

rule basecalling:
    input:
        'flags/0_prepare.done'
    output: touch("flags/1_basecalling.done")
    log:
        log = 'log/1_basecalling.log'
    params:
        basecaller = config['basecaller'],
        cfg = config,
        cmdline = 'log/cmdline.log',
    run:
        if params.basecaller == "dorado" and config["info_dict"]['barcoding']:
            print('[red] dorado is not yet compatible with barcoding [/red]')
            sys.exit(1)

        while not gpu_available():
            time.sleep(60)

        if params.basecaller == "dorado":
            print('[yellow] run dorado [/yellow]')
            cmd = dorado_basecalling(
                params.cfg,
                params.cmdline,
                log.log
            )
        elif params.basecaller == "guppy":
            print('[yellow] run guppy [/yellow]')  
            cmd = guppy_basecalling(
                params.cfg,
                params.cmdline,
                log.log
            )
        else:
            print('[red]basecaller unknown: {} [/red]'.format(params.basecaller))
            sys.exit(1)
