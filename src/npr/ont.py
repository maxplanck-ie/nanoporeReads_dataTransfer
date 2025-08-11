import os
import signal
import sys
from importlib.metadata import version
from pathlib import Path
import subprocess as sp
from threading import Event
import rich_click as click
import yaml
from rich import print
from npr.communication import query_parkour, send_email, ship_qcreports, standard_text
from npr.ont_pipeline import find_new_flowcell, get_periphery, read_flowcell_info, read_samplesheet
from npr.snakehelper import getfast5foot, get_disk_stat, scan_multiqc, config_to_smkcmd, print_header


# set up CLI args.
@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-c",
    "--configfile",
    type=click.Path(exists=True),
    required=False,
    default=os.path.expanduser("~/configs/ont_prod.yaml"),
    help="Specify a custom yaml file.",
    show_default=True,
)
@click.option(
    "-d",
    "--directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=False,
    default=None,  # Default: unspecified (defined in config)
    help="Specify a custom input directory.",
)
@click.option(
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Run snakemake as --dryrun",
)
@click.option(
    "--organism",
    default=None,
    show_default=True,
    help="If known, specify organism as --organism",
)
@click.option(
    "--protocol",
    default=None,
    show_default=True,
    help="If known, specify protocol as --protocol",
)
@click.option(
    "-f",
    "--flowcell",
    default=None,
    show_default=True,
    help="Target a specific flowcell. Can be a substring of the directory to look for in the offload directory (config, or via --directory).",
)
@click.option(
    "--force",
    default=False,
    is_flag=True,
    show_default=True,
    help="If a flowcell is specified, force it to run (even if the appropriate flag is not set in that folder).",
)

# run workflow.
def ont(**kwargs):
    # Crash of force and no flowcell.
    if kwargs["force"] and not kwargs["flowcell"]:
        print("Force flag set without specifying a flowcell. Exiting.")
        sys.exit(1)

    # print what config is used.
    print(
        "Starting pipeline with config: [green]{}[/green]".format(kwargs["configfile"])
    )

    # Load config from file
    config = yaml.safe_load(open(kwargs["configfile"]))

    # update config if runtime args have been set
    if kwargs["directory"] is not None:
        config["paths"]["offloadDir"] = kwargs["directory"]
    print("Watching directory: [green]{}[/green]".format(config["paths"]["offloadDir"]))

    if kwargs["dryrun"] is not False:
        config["snakemake"]["dryrun"] = True
        if config["options"]["verbosity"]:
            print("Use snakemake in dryrun.")

    if kwargs["flowcell"] is not False:
        config["target_flowcell"] = kwargs["flowcell"]
        config["force_processing"] = kwargs["force"]
        if config["options"]["verbosity"]:
            print(
                "Target flowcell is [green]{}[/green].".format(
                    config["target_flowcell"]
                )
            )

    # snakefile to config['snakemake']
    config["snakemake"]["snakefile"] = os.path.join(__file__, 'rules_dorado', 'ont_pipeline.smk')

    # initialize config['info_dict']
    # this applies only to the _first_ flowcell (used to sidetrack Parkour query)
    # notice that info_dict is flowcell-specific and will be reset/re-evalutated for each
    if "info_dict" not in config:
        config["info_dict"] = {}
    if kwargs["organism"] is not None:
        config["info_dict"]["organism"] = kwargs["organism"]
        print("Set organism to {}".format(config["info_dict"]["organism"]))
    if kwargs["protocol"] is not None:
        config["info_dict"]["protocol"] = kwargs["protocol"]
        print("Set protocol to {}".format(config["info_dict"]["protocol"]))
    # ensure that pipeline version is tracked in metadata.yml
    config["info_dict"]["pipeline_version"] = version("npr")

    # add conda env path where executable live
    config["paths"]["conda_env"] = sys.prefix

    # start workflow.
    main(config)


def main(config):
    while True:
        # Set HUP
        HUP = Event()

        def breakSleep(signo, _frame):
            HUP.set()

        def sleep():
            HUP.wait(timeout=float(60 * 60 * int(config["options"]["sleep_time"])))

        signal.signal(signal.SIGHUP, breakSleep)

        flowcell, msg, base_path = find_new_flowcell(config)
        if flowcell:
            print_header(flowcell)
            print_header('parkour')
            msg = query_parkour(config, flowcell, msg)
            # The following should be simplified but I did not want to touch
            # read_flow_cell_info() for nowq
            print_header('read_flowcell_info')
            config["info_dict"] = read_flowcell_info(
                config, config["info_dict"], base_path
            )

            # read samplesheet
            print_header('read_samplesheet')
            bc_kit, data = read_samplesheet(config)
            config["data"] = data
            config["bc_kit"] = bc_kit
            print(config["data"])
            print("samplesheet is read sucessfully")

            # Snakemake setup.
            print_header('Starting snakemake')
            config["info_dict"]["logfile"] = os.path.join(
                config["info_dict"]["flowcell_path"], "log", "ont.log"
            )

            # write the updated config file under the output path
            configFile = os.path.join(
                config["info_dict"]["flowcell_path"], "pipeline_config.yaml"
            )
            config["info_dict"]["configFile"] = configFile
            with open(configFile, "w") as f:
                yaml.dump(config, f, default_flow_style=False)

            config["info_dict"]["transfer_path"] = get_periphery(config)

            print("[green]Starting snakemake. [/green]")
            print("file    {}".format(config["snakemake"]["snakefile"]))
            print("config  {}".format(config["info_dict"]["configFile"]))
            print("workdir {}".format(config["info_dict"]["flowcell_path"]))
            print("transdir {}".format(config["info_dict"]["transfer_path"]))

            # static parameters are defined in dict config['snakemake']
            # flowcell specific parameters are taken from config['info_dict']
            _smk_cmd = config_to_smkcmd(config["snakemake"])
            _smk_cmd.extend(
                [
                    '--configfile', config["info_dict"]["configFile"],
                    '-d', config["info_dict"]["flowcell_path"],
                ]
            )
            _smk_cmd.insert(0, 'snakemake')
            
            #
            print(f"Starting snakemake with command:\n {' '.join(_smk_cmd)}\n")
            rg = sp.Popen(_smk_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, text=True, bufsize=1, universal_newlines=True)
            ols = []
            for line in iter(rg.stdout.readline, ''):
                print(line, end='')
                ols.append(line)
            rg.stdout.close()
            smk_log = ''.join(ols)
            return_code = rg.wait()

            #stdout, stderr = rg.communicate()
            if return_code != 0:
                msg += f"Snakemake errorcode {return_code}\n"
                msg += smk_log
                send_email("Snakemake failed for flowcell:", msg, config, failure=True)
                sys.exit(1)
            else:
                with open(os.path.join(config["info_dict"]["flowcell_path"], "snakemake.log"), 'w') as f:
                    f.write(smk_log)
            print(config)
            print_header("Post snakemake processing")

            if return_code == 0 and not config["snakemake"]["dryrun"]:
                msg += "pod5 compression: {}\n".format(
                    getfast5foot(
                        config["info_dict"]["base_path"],
                        config["info_dict"]["flowcell_path"],
                    )
                )
                # spread the news
                ship_qcreports(config, flowcell)
                config["QC"] = scan_multiqc(config)
                config["SM"] = {
                                    "offload (deepseq)": get_disk_stat(config['paths']['offloadDir']),
                                    "inbox (bioinfo)": get_disk_stat(config['paths']['outputDir']),
                                    "periphery": get_disk_stat(config['info_dict']['transfer_path'], True)
                                }
            
                msg = standard_text(config)
                send_email("Successfully finished flowcell:", msg, config, failure=False)

                summaryFile = os.path.join(
                    config["info_dict"]["flowcell_path"], "log/summary.yaml"
                )
                with open(summaryFile, "w") as f:
                    yaml.dump(config, f, default_flow_style=False)

                Path(
                    os.path.join(config["info_dict"]["flowcell_path"], "analysis.done")
                ).touch()

                # wipe config['info_dict']
                config["info_dict"] = {}
                
        else:
            print("No flowcells found. I go back to sleep.")
            sleep()
            continue
