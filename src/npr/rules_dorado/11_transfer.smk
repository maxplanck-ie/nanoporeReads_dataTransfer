import time
import yaml
import shutil
from npr.ont_pipeline import get_dest_path

rule transfer_11:
    input:
        expand("transfer/Project_{sample_project}/QC/multiqc_report.html", sample_project=set(sample_projects))
    output:
        touch("transfer.done"),
    log:
        file="log/11_transfer.log"
    run:
        dirs = glob.glob(os.path.join("transfer", 'Project*'))
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
                os.chmod(dest, 0o700)
                
                # write metadata
                metayaml = os.path.join(dest, 'metadata.yaml') 
                with open(metayaml, 'w') as yaml_file:
                    yaml.dump(config['info_dict'], yaml_file, default_flow_style=False)

                # Include a transfer of the 'raw reports'
                reports_dir = os.path.join(dest, "raw_sequencer_reports")
                if not os.path.exists(reports_dir):
                    os.makedirs(reports_dir)
                    os.chmod(reports_dir, 0o700)
                for raw_report in os.listdir("reports"):
                    inf = os.path.join("reports", raw_report)
                    outf = os.path.join(reports_dir, raw_report)
                    if os.path.isfile(inf):
                        shutil.copy2(inf, outf)
                        os.chmod(outf, 0o700)

                # strip user write permission from destination directory 'dest'
                for r, dirs, files in os.walk(dest):
                    for d in dirs:
                        os.chmod(os.path.join(r, d), 0o700)
                    for f in files:
                        os.chmod(os.path.join(r, f), 0o700)
                log.write(f"{ts}: Stripped user write permission from {dest}\n")