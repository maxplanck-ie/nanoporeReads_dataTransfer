from npr.snakehelper import fast5_to_pod5
from npr.snakehelper import merge_pod5
from npr.snakehelper import getfast5foot

rule prepare:
    output:
       touch('flags/0_prepare.done')
    params:
        idir = config["info_dict"]["base_path"],
        baseout = config['info_dict']['flowcell_path'],
        log  = 'log/cmdline.log'
    run:
        if os.path.exists(os.path.join(params.baseout, "pod5", "merged.pod5")):
            print("Merged pod5 file exist")
        elif os.path.exists(os.path.join(params.idir, "pod5_pass")):
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

