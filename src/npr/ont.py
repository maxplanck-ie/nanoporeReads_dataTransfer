# CLI / pretty print.
from rich import print, inspect
import rich_click as click
# default libs
import glob
import os
import snakemake
import sys
import yaml
# sleeper libs
from threading import Event
import signal
# npr libs
from npr.ont_pipeline import find_new_flowcell
from npr.ont_pipeline import read_flowcell_info
from npr.ont_pipeline import read_samplesheet
from npr.ont_pipeline import config_to_basecallcmd
from npr.communication import query_parkour, send_email
import subprocess as sp
from importlib.metadata import version
from pathlib import Path


# set up CLI args.
@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"]
    )
)
@click.option(
   "-c",
   "--configfile",
   type=click.Path(exists=True),
   required=False,
   default=os.path.expanduser('~/configs/ont_prod.yaml'),
   help='specify a custom yaml file.',
   show_default=True
)
# run workflow.

def ont(configfile):
    # print what config is used.
    print(
        "Starting pipeline with: [green]{}[/green]".format(configfile)
    )
    # Load config up.
    config = yaml.safe_load(
        open(configfile)
    )
    # start workflow.
    main(config)

def main(config):
    while True:
        # Set HUP
        HUP = Event()
        def breakSleep(signo, _frame):
            HUP.set()
        def sleep():
            HUP.wait(timeout=float(60*60))
        signal.signal(signal.SIGHUP, breakSleep)

        flowcell, msg = find_new_flowcell(config)
        if flowcell:
            info_dict, msg = query_parkour(config, flowcell, msg)
            config["info_dict"] = read_flowcell_info(config, info_dict)

            # read samplesheet
            bc_kit,data = read_samplesheet(config)
            config["data"] = data
            config["bc_kit"] = bc_kit
            print(config["data"])
            print("samplesheet is read sucessfully")
            # Include the rulePath in the configFile.
            config['paths']['rulesPath'] = os.path.join(
                os.path.realpath(os.path.dirname(__file__)),
                'rules'
            )
            # Infer the basecallcmd already.
            
            config['basecallcmd'] = config_to_basecallcmd(config)
            # write the updated config file under the output path
            configFile = os.path.join(
                config["paths"]["outputDir"],
                config["input"]["name"],
                "pipeline_config.yaml"
            )
            with open(configFile, 'w') as f:
                yaml.dump(config, f, default_flow_style=False)

            #run snakemake
            output_directory = os.path.join(
                config["paths"]["outputDir"],
                config["input"]["name"]
            )
            snakefile_file = os.path.join(
                os.path.realpath(os.path.dirname(__file__)),
                "rules",
                "ont_pipeline.smk"
            )

            # snakemake log file
            fnames = glob.glob(
                os.path.join(output_directory, 'ont_run-[0-9]*.log')
            )
            if len(fnames) == 0:
                n = 1  # no matching files, then this is the first run
            else:
                fnames.sort(key=os.path.getctime)
                n = int(fnames[-1].split("-")[-1].split(".")[0]) + 1  # get new run number
            
            snak_stat = snakemake.snakemake(
                snakefile = snakefile_file,
                debug = True,
                #cores = 25,
                max_jobs_per_second = 1,
                printshellcmds = True,
                verbose = True,
                configfiles = [configFile],
                workdir = output_directory
            )
            if not snak_stat:
                msg += "snake crashed."
                sys.exit()
            Path(os.path.join(output_directory, 'analysis.done')).touch()
            print(config)
            msg = 'guppy version: {}\n'.format(config['guppy_basecaller']['guppy_version'])
            msg += 'guppy model: {}\n'.format(config['guppy_model_used'].split('/')[-1])
            msg += 'minimap2 version: {}\n\n'.format(config['mapping']['minimap2_version'])
            msg += "flowcell {} is analysed successfully".format(flowcell)
            send_email(msg, version('npr'), os.path.basename(flowcell), config)
        else:
            print("No flowcells found. I go back to sleep.")
            sleep()
            continue