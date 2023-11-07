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
from npr.ont_pipeline import get_periphery
from npr.communication import query_parkour, send_email, ship_qcreports, standard_text
from npr.snakehelper import getfast5foot, get_seqdir, scan_multiqc, monitor_storage
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

@click.option(
    "-d",
    "--directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=False,
    default=None,  # Default: unspecified (defined in config)
    help="Specify a custom input directory."
)

@click.option(
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Run snakemake as --dryrun"
)

@click.option(
    "--organism",
    default=None,
    show_default=True,
    help="If known, specify organism as --organism"
)

@click.option(
    "--protocol",
    default=None,
    show_default=True,
    help="If known, specify protocol as --protocol"
)
# run workflow.
def ont(**kwargs):
    # print what config is used.
    print(
        "Starting pipeline with config: [green]{}[/green]".format(kwargs['configfile'])
    )

    # Load config from file
    config = yaml.safe_load(open(kwargs['configfile']))

    # update config if runtime args have been set
    if ( kwargs['directory'] is not None):
        config['paths']['offloadDir'] = kwargs['directory']
    print(
        "Watching directory: [green]{}[/green]".format(config['paths']['offloadDir'])
    )

    if ( kwargs['dryrun'] is not False):
        config['snakemake']['dryrun'] = True
        print("Use snakemake in dryrun.")

    # add rulesPath to config['paths'] _not_ to config['snakemake']
    # since 'rulesPath' is not a snakemake option
    config['paths']['rulesPath'] = os.path.join(
        os.path.realpath(os.path.dirname(__file__)),
        config['paths']['rulesDir']
    )

    # snakefile to config['snakemake']
    config['snakemake']['snakefile'] = os.path.join(
        config['paths']['rulesPath'],
        "ont_pipeline.smk"
    )

    # initialize config['info_dict']
    # this applies only to the _first_ flowcell (used to sidetrack Parkour query)
    # notice that info_dict is flowcell-specific and will be reset/re-evalutated for each
    if 'info_dict' not in config:
        config['info_dict'] = {}
    if ( kwargs['organism'] is not None):
        config['info_dict']['organism'] = kwargs['organism']
        print("Set organism to {}".format(config['info_dict']['organism']))
    if ( kwargs['protocol'] is not None):
        config['info_dict']['protocol'] = kwargs['protocol']
        print("Set protocol to {}".format(config['info_dict']['protocol']))
    # add conda env path where executable live
    config['paths']['conda_env'] = sys.prefix

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


        flowcell, msg, base_path = find_new_flowcell(config)
        if flowcell:
            if ('organism' not in config['info_dict'] or
                'protocol' not in config['info_dict']):
                # need parkour query only if 'organism' or 'protocol' is undefined
                msg = query_parkour(config, flowcell, msg)
            # The following should be simplified but I did not want to touch
            # read_flow_cell_info() for nowq
            config["info_dict"] = read_flowcell_info(config, config["info_dict"], base_path)

            # read samplesheet
            bc_kit,data = read_samplesheet(config)
            config["data"] = data
            config["bc_kit"] = bc_kit
            print(config["data"])
            print("samplesheet is read sucessfully")

            config['info_dict']['logfile']=os.path.join(
                config['info_dict']['flowcell_path'],
                "log",
                "ont.log"
            )

            # write the updated config file under the output path
            configFile = os.path.join(
                config['info_dict']['flowcell_path'],
                "pipeline_config.yaml"
            )
            config['info_dict']['configFile']=configFile
            #print(config)
            with open(configFile, 'w') as f:
                yaml.dump(config, f, default_flow_style=False)

            send_email("Found flowcell:", msg, config, allreceivers=False)

            """
            # This seems to be unused
            output_directory = config['info_dict']['flowcell_path']
            # snakemake log file
            fnames = glob.glob(
                os.path.join(output_directory, 'ont_run-[0-9]*.log')
            )
            if len(fnames) == 0:
                n = 1  # no matching files, then this is the first run
            else:
                fnames.sort(key=os.path.getctime)
                n = int(fnames[-1].split("-")[-1].split(".")[0]) + 1  # get new run number
            """
            config['info_dict']['transfer_path'] = get_periphery(config)

            print("[green]Starting snakemake. [/green]")
            print("file    {}".format(config['snakemake']['snakefile']))
            print("config  {}".format(config['info_dict']['configFile']))
            print("workdir {}".format(config['info_dict']['flowcell_path']))
            print("transdir {}".format(config['info_dict']['transfer_path']))

            # static parameters are defined in dict config['snakemake']
            # flowcell specific parameters are taken from config['info_dict']
            snak_stat = snakemake.snakemake(
                **config['snakemake'],
                configfiles = [ config['info_dict']['configFile'] ],
                workdir = config['info_dict']['flowcell_path'],
            )
            if not snak_stat:
                msg += "snake crashed with {}".format(snak_stat)
                send_email("Snakemake failed for flowcell:", msg, config)
                sys.exit(1)

            msg += 'pod5 compression: {}\n'.format(
                getfast5foot(
                    config['info_dict']['base_path'],
                    config['info_dict']['flowcell_path']
                )
            )
            # spread the news
            ship_qcreports(config, flowcell)
            config['QC'] = scan_multiqc(config)
            config['SM'] = monitor_storage(config)
            msg = standard_text(config)
            send_email("Successfully finished flowcell:", msg, config)

            summaryFile = os.path.join(
                config['info_dict']['flowcell_path'],
                "log/summary.yaml"
            )
            with open(summaryFile, 'w') as f:
                yaml.dump(config, f,default_flow_style=False)

            Path(os.path.join(config['info_dict']['flowcell_path'], 'analysis.done')).touch()

            # wipe config['info_dict']
            config['info_dict']={}
        else:
            print("No flowcells found. I go back to sleep.")
            sleep()
            continue

