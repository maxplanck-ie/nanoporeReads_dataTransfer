from npr.snakehelper import basecalling

rule guppy_basecalling:
    input:
        'flags/0_prepare.done'
    output: touch("flags/1_basecalling.done")
    log:
        log = 'log/basecalling.log'
    params:
        cfg = config,
        cmdline = 'log/cmdline.log'
    run:
        cmd = basecalling(
            params.cfg,
            params.cmdline,
            log.log
        )
