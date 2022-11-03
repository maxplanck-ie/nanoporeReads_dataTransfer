from npr.snakehelper import basecalling
from npr.snakehelper import fast5_to_pod5

rule pod5:
    output:
       touch('flags/1_pod5.done')
    params:
        idir = config["info_dict"]["base_path"],
        baseout = config['info_dict']['flowcell_path']
    run:
        fast5_to_pod5(
            params['idir'],
            params['baseout']
        )

rule basecall:
    input:
        'flags/1_pod5.done'
    output: touch("flags/1_basecalling.done")
    params:
        cfg = config
    run:
        basecalling(params['cfg'])