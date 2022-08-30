
rule basecall:
    output: touch("demux.done")
    params:
        basecallcmd = config['basecallcmd']
    log:
        out = "LOG/guppy.log.out",
        err = "LOG/guppy.log.err"
    shell:'''
    {params.basecallcmd} > {log.out} 2> {log.err}
    '''