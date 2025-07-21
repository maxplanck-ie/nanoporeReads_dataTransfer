import time
import yaml
from npr.ont_pipeline import get_dest_path

rule transfer:
    input:
        flag_mulitqc="flags/08_multiqc.done",
        flag_modbed=lambda wildcards: "flags/07_modbed.done" if do_modbed and do_align else []
    output:
        touch("flags/09_transfer.done"),
    log:
        file="log/09_transfer.log"
    run:
        dirs = glob.glob(os.path.join(transfer_dir, 'Project*'))
        with open(log.file, "a") as log:
            if not dirs:
                ts = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                log.write(f"{ts}: No directory found.\n")
            for dir in dirs:
                # timestamp
                ts = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                if not os.path.isdir(dir):
                    log.write(f"{ts}: Do not transfer {dir}\n")
                    continue
                # get destination
                dest= get_dest_path(config,dir)
                log.write(f"{ts}: mkdir {dest}; cp -r {dir} {dest}\n")
                # copy data over
                shell("mkdir -p {dest}; cp -r {dir} {dest}")
                # write metadata
                metayaml = os.path.join(dest, 'metadata.yaml') 
                with open(metayaml, 'w') as yaml_file:
                    yaml.dump(config['info_dict'], yaml_file, default_flow_style=False)

                # strip user write permission from destination directory 'dest'
                for r, dirs, files in os.walk(dest):
                    for d in dirs:
                        os.chmod(os.path.join(r, d), 0o700)
                    for f in files:
                        os.chmod(os.path.join(r, f), 0o700)

                log.write(f"{ts}: Stripped user write permission from {dest}\n")
