from npr.snakehelper import basecalling
from npr.snakehelper import fast5_to_pod5
from npr.snakehelper import merge_pod5

rule prepare_pod5:
    output:
       touch('flags/1_pod5.done')
    params:
        idir = config["info_dict"]["base_path"],
        baseout = config['info_dict']['flowcell_path'],
        log  = 'log/cmdline.log'
    run:
        if os.path.exists(os.path.join(params.idir, "pod5_pass")):
            merge_pod5(
                params.idir, 
                params.baseout, 
                params.log)
        elif os.path.exists(os.path.join(params.idir, "fast5_pass")):
            fast5_to_pod5(
                params.idir,
                params.baseout,
                params.log
            )

rule guppy_basecalling:
    input:
        'flags/1_pod5.done'
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
