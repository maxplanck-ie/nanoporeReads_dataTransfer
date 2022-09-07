from npr.snakehelper import config_to_basecallcmd

rule basecall:
    output: touch("flags/1_basecalling.done")
    params:
        basecallcmd = config_to_basecallcmd(config)
    log:
        out = "log/guppy.log.out",
        err = "log/guppy.log.err"
    shell:'''
    {params.basecallcmd} > {log.out} 2> {log.err}
    '''
