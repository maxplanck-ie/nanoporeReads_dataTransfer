# ont.py
# Main script to monitor new data and launch new jobs
# (C) 2024 Bioinformatics Core
# Max Plank Institute for Immunobiology and Epigenetics

import glob

# default libs
import os

# sleeper libs
import signal
import sys
from importlib.metadata import version
from pathlib import Path
from threading import Event

#libs for asynchronous functions
import asyncio, asyncssh

import rich_click as click
import snakemake
import yaml

# CLI / pretty print.
from rich import print

from npr.communication import query_parkour, send_email, ship_qcreports, standard_text, transfer_to_remote

# npr libs
from npr.ont_pipeline import (
    find_new_flowcell,
    get_periphery,
    read_flowcell_info,
    read_samplesheet,
    sanitize_info_dict_for_remote
)
from npr.snakehelper import getfast5foot, monitor_storage, scan_multiqc, merge_dicts


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
    "--send_to_remote",
    is_flag=True,
    default=False,
    show_default=True,
    help="Send pod5 files for a specific flowcell to a remote vm. Requires flowcell, organism and protocol to be specified.",
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
    #If send_to_remote specified, require flowcell and organism. Else exit.
    if kwargs["send_to_remote"] and not all([kwargs["flowcell"],kwargs["organism"],kwargs["protocol"]]):
        print("Sending to remote requires flowcell, organism and protocol specification.")
        sys.exit(1)

    # print what config is used.
    print(
        "Starting pipeline with config: [green]{}[/green]".format(kwargs["configfile"])
    )

    # Load config from file
    config = yaml.safe_load(open(kwargs["configfile"]))

    #If running on remote, exit if send_to_remote is requested
    if kwargs["send_to_remote"] and config["remote_vm"]["is_remote"]:
        print("The workflow is running on a remote machine. Sending to remote is not possible.")
        sys.exit(1)

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

    # add rulesPath to config['paths'] _not_ to config['snakemake']
    # since 'rulesPath' is not a snakemake option
    config["paths"]["rulesPath"] = os.path.join(
        os.path.realpath(os.path.dirname(__file__)), config["paths"]["rulesDir"]
    )

    # snakefile to config['snakemake']
    config["snakemake"]["snakefile"] = os.path.join(
        config["paths"]["rulesPath"], "ont_pipeline.smk"
    )

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
    
    config["options"]["send_to_remote"] = kwargs["send_to_remote"]

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
            if not config["remote_vm"]["is_remote"]:
                if (
                    "organism" not in config["info_dict"]
                    or "protocol" not in config["info_dict"]
                ):
                    # need parkour query only if 'organism' or 'protocol' is undefined
                    msg = query_parkour(config, flowcell, msg)
                # The following should be simplified but I did not want to touch
                # read_flow_cell_info() for nowq
                config["info_dict"] = read_flowcell_info(
                    config, config["info_dict"], base_path
                )

                # read samplesheet
                bc_kit, data = read_samplesheet(config)
                config["data"] = data
                config["bc_kit"] = bc_kit
                print(config["data"])
                print("samplesheet is read sucessfully")

            else:
               #update current config with remote
               remote_pipeline_config_file = os.path.join(base_path,"remote_pipeline_config.yaml")
               remote_pipeline_config = yaml.safe_load(open(remote_pipeline_config_file))
               #config=merge_dicts(remote_pipeline_config,config) 
               config = remote_pipeline_config
               config = sanitize_info_dict_for_remote(config,flowcell,basepath)
            #same for local and remote
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

            send_email("Found flowcell:", msg, config, allreceivers=False)

            if config["options"]["send_to_remote"]:
               # Path(
               #     os.path.join(config["info_dict"]["flowcell_path"], "ignore.on.local")
               # ).touch()
                remote_config=yaml.safe_load(open(config["remote_vm"]["remote_config"]))
                remote_pipeline_config=merge_dicts(config,remote_config)
                remote_configFile = os.path.join(config["info_dict"]["flowcell_path"], "remote_pipeline_config.yaml")
                with open(remote_configFile, "w") as f:
                    yaml.dump(remote_pipeline_config, f, default_flow_style=False)
                try:
                    msg = asyncio.run(transfer_to_remote(flowcell,config,msg))
                except (OSError, asyncssh.Error) as exc:
                    sys.exit('SSH transfer failed: ' + str(exc))
                msg += "Transfer to remote completed successfully. \n"
                send_email("Successfully finished flowcell:", msg, config)
                exit(0)

            else:
                config["info_dict"]["transfer_path"] = get_periphery(config)

                print("[green]Starting snakemake. [/green]")
                print("file    {}".format(config["snakemake"]["snakefile"]))
                print("config  {}".format(config["info_dict"]["configFile"]))
                print("workdir {}".format(config["info_dict"]["flowcell_path"]))
                print("transdir {}".format(config["info_dict"]["transfer_path"]))

                # static parameters are defined in dict config['snakemake']
                # flowcell specific parameters are taken from config['info_dict']
                snak_stat = snakemake.snakemake(
                    **config["snakemake"],
                    configfiles=[config["info_dict"]["configFile"]],
                    workdir=config["info_dict"]["flowcell_path"],
                )
                if not snak_stat:
                    msg += f"snake crashed with {snak_stat}"
                    send_email("Snakemake failed for flowcell:", msg, config)
                    sys.exit(1)

                msg += "pod5 compression: {}\n".format(
                    getfast5foot(
                        config["info_dict"]["base_path"],
                        config["info_dict"]["flowcell_path"],
                    )
                )
            if not config["remote_vm"]["is_remote"]:
                # spread the news
                ship_qcreports(config, flowcell)
                config["QC"] = scan_multiqc(config)
                config["SM"] = monitor_storage(config)
                msg = standard_text(config)
                send_email("Successfully finished flowcell:", msg, config)

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
