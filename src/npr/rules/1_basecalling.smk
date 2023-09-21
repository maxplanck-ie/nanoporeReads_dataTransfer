from npr.snakehelper import guppy_basecalling, dorado_basecalling
from rich import print

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
            params.basecaller = "guppy"
            print('[red] dorado is not yet compatible with barcoding [/red]')
            print('[red] Overwrite params.basecaller: {} [/red]'.format(params.basecaller))

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

""" 
rule guppy_basecalling:
    input:
        'flags/0_prepare.done'
    output: touch("flags/1_basecalling.done")
    log:
        log = 'log/1_basecalling.log'
    params:
        cfg = config,
        cmdline = 'log/cmdline.log'
    run:
        cmd = basecalling(
            params.cfg,
            params.cmdline,
            log.log
        )
 """