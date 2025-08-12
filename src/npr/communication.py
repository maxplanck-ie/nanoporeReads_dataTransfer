import glob
import os
import shutil
import smtplib
import sys
import json
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from importlib.metadata import version
from time import sleep
import requests
from rich import print
from pathlib import Path
from dominate.tags import html, div, br, b
from tabulate import tabulate
from rich import print

def ship_qcreports(config, flowcell):
    """
    Legacy shipment of QC reports. Reconsider since this generates quite some duplication
    The idea is to solve
     1. NGS team can only reads via SAMBA share
     2. bioinfo may not be able to open html on analysis server
    """
    # Create flowcell folder if it doesn't exist.
    yrstr = flowcell[:4]
    sambadir = Path(config['paths']['deepseq_qc'], f"Sequence_Quality_{yrstr}", f"ONT_{yrstr}")
    sambadir.mkdir(parents=True, exist_ok=True)

    fc_dir = Path(config["info_dict"]["flowcell_path"])

    samba_target = sambadir / fc_dir.name
    bioinfo_target =Path(config["paths"]["bioinfocoredir"], fc_dir.name)
    print("[green]Copy QC reports[/green]")
    print(f"... to bioinfo: {bioinfo_target}")
    print(f"... to sambahost: {samba_target}")

    # copy reports
    reports = glob.glob(
        os.path.join(config["info_dict"]["base_path"], "reports", "*.html")
    )
    for report in reports:
        print(f"copying... {report}")
        shutil.copy(report, bioinfo_target)
        shutil.copy(report, samba_target)

    # copy QC directories: transfer/Project*/QC
    for local_dir in glob.glob(f"{fc_dir}/transfer/Project*/QC/"):
        p_dir = os.path.dirname(os.path.dirname(local_dir))
        p_name = os.path.basename(p_dir)
        target_dir = os.path.join(fc_dir.name, p_name, "QC")

        # copy to bioinfo
        bioinfo_target = os.path.join(config["paths"]["bioinfocoredir"], target_dir)
        print(f"copy {local_dir}")
        shutil.copytree(local_dir, bioinfo_target, dirs_exist_ok=True)
        shutil.copytree(local_dir, samba_target, dirs_exist_ok=True)


def standard_text(config):
    """
    Draft a short letter to end user that will be part of the success email
    Include also QC metrics obtained from multiqc report
    """

    # Relevant keys from info_dict for email
    relkeys = ['pipeline_version', 'model', 'parkour_protocol', 'modifications', 'organism', 'flowcell', 'kit', 'barcoding', 'barcode_kit']
    # Get the project(s)
    pid_to_fids = query_parkour_project(config)

    # Create the standard succes email
    msg = f'Project: {config["data"]["projects"]}\n'
    for pid in pid_to_fids:
        msg += f'Expected flowcells for {pid}: {len(pid_to_fids[pid])}\n'
    msg += f'Output Folder: {config["info_dict"]["transfer_path"]}\n' 
    msg += f'Protocol: {config["info_dict"]["parkour_protocol"]}\n'
    msg += f'Flowcell: {config["info_dict"]["flowcell"]}\n'
    msg += f'Kit: {config["info_dict"]["kit"]}\n'
    msg += f'Barcoding: {config["info_dict"]["barcoding"]}\n'
    msg += f'Barcode Kit: {config["info_dict"]["barcode_kit"]}\n'
    # Push in the storage settings:
    for disk in config["SM"]:
        msg += f'Storage occupied {disk}: {config["SM"][disk]["percentage"]}%\n'

    return msg


def send_email(subject, body, config, failure=False):
    """
    Send email including key information about the run
    Also print message to stdout
    """
    # Set up mailer.
    mailer = MIMEMultipart("alternative")
    _fid = os.path.basename(config['info_dict']['base_path'])
    mailer["Subject"] = f"[npr] [{version("npr")}] {subject} {_fid}"    
    mailer["From"] = config["email"]["from"]
    if failure:
        _receivers = config["email"]["failure"].split(',')
    else:
        _receivers = config["email"]["to"].split(',')
    mailer["To"] = ", ".join(_receivers)

    # body comes in as a string (\n delimited.)
    # Convert to HTML to be pretty.
    _html = html()
    for _l in body.splitlines():
        if ':' in _l:
            label, value = _l.split(":", 1)
            _html.add(div(b(label + ":"), " " + value, br()))
        else:
            _html.add(div(_l, br()))
    _html.add(br())
    if 'QC' in config:
        if 'GI' in config['QC']:
            for metric in config['QC']['GI']:
                _html.add(div(b(metric + ":"), " " + f"{config['QC']['GI'][metric]}", br()))
            _html.add(br())
    
    if 'QC' in config and 'QC' in config['QC']:
        _tablehead = ['Sample name', 'Project', 'Kraken top hit', 'Parkour organism', 'Total bp', 'Total reads', 'N50', 'Median length', 'Median Q', '% reads Q >= 18']
        _tablecont = []
        for sid in config['QC']['QC']:
            _tablecont.append(
                [sid] + [config['QC']['QC'][sid][k] for k in _tablehead]
            )
        _html = _html.render() + tabulate(_tablecont, ['Sample ID'] + _tablehead, tablefmt="html")
    else:
        _html = _html.render()

    email = MIMEText(_html, 'html')
    mailer.attach(email)
    if config["email"]["host"] is not None:
        s = smtplib.SMTP(config["email"]["host"])
        s.sendmail(config["email"]["from"], _receivers, mailer.as_string())


def query_parkour(config, flowcell, msg):
    """
    query parkour to update info_dict (protocol and organism)
    """

    if not config["parkour"]["url"]:
        msg += "Parkour URL not specified."
        return msg

    try:
        fc = flowcell.split("_")[3]
    except Exception as e:
        print(f"[red]Exception: {e}[/red]")
        msg += 'Parkour Error: flowcell "{}" cannot be queried in parkour.'.format(
            flowcell
        )
        send_email("Error with flowcell:", msg, config, failure=False)
        sys.exit(1)
    # test for flow cell re-use.
    # flowcell that's re-used gets higher increment.
    postfixes = ["-5", "-4", "-3", "-2", "-1", ""]
    flowcellqueries = []
    for pf in postfixes:
        d = {"flowcell_id": fc + pf}
        flowcellqueries.append(fc + pf)
        res = requests.get(
            config["parkour"]["url"] + "/api/analysis_list/analysis_list/",
            auth=(config["parkour"]["user"], config["parkour"]["password"]),
            params=d,
            verify=config["parkour"]["pem"],
        )
        if res.status_code == 200:
            msg += f"Parkour query 200 for fid {fc} with re-use index {pf}.\n"
            msg += "\n"
            msg += "Parkour queried with:\n"
            for fq in flowcellqueries:
                msg += f"{fq}\n"
            msg += "\n"
            parkour_dict = res.json()
            print(parkour_dict)
            first_key = list(parkour_dict.keys())[0]
            first_entry = list(parkour_dict[first_key].keys())[0]
            # Parkour pulls give back a list containing [name, label, yaml]
            # In the config the name: [fna, bed] is encoded.
            organism_tupe = parkour_dict[first_key][first_entry][3]
            protocol = parkour_dict[first_key][first_entry][2]

            if protocol not in config['dorado_basecaller']['protocol_models']:
                print(f"[red]Protocol {protocol} not found in config. Setting default {config['dorado_basecaller']['default_model']}[/red]")
                config["info_dict"]["model"] = config['dorado_basecaller']['default_model']
            else:
                print(f"Setting model to {config['dorado_basecaller']['protocol_models'][protocol]}")
                config["info_dict"]["model"] = config['dorado_basecaller']['protocol_models'][protocol]
            config['info_dict']['parkour_protocol'] = protocol
            # Note that at this stage it's simple to extract if we have modification calls or not:
            # Modifications will have a ',' in the 'model' field.
            config["info_dict"]["modifications"] = ("," in config["info_dict"]["model"])
            
            if (
                fc == "PAK78871"
                or fc == "PAK79330"
                or fc == "PAK77043"
                or fc == "FAV22714-2"
                or fc == "PAK78965"
                or fc == "PAK79757"
            ):
                print("[red] Protocol override [/red]")
                protocol = "cdna"
            if "cDNA" in protocol:
                protocol = "cdna"
            elif "DNA" in protocol:
                protocol = "dna"
            elif "RNA" in protocol:
                protocol = "rna"
            else:
                print(protocol)
                print("protocol not found. Default to dna.")
                protocol = "dna"

            # update config['info_dict']
            config["info_dict"]["protocol"] = protocol
            
            print(f"Parkour organism received = {organism_tupe}")

            if organism_tupe[0] not in config["genome"].keys():
                config["info_dict"]["organism"] = "other"
                config["info_dict"]["organism_genome"] = None
                config["info_dict"]["organism_genes"] = None
                config["info_dict"]["organism_label"] = None
            else:
                _org = organism_tupe[0]
                config["info_dict"]["organism"] = _org
                # encoded in config as list of [genome.fa, genes.bed]
                config["info_dict"]["organism_genome"] = config["genome"][_org][0]
                config["info_dict"]["organism_genes"] = config["genome"][_org][1]
                config["info_dict"]["organism_label"] = organism_tupe[1]
            return msg

    msg += f"Parkour query failed for {fc}.\n"
    msg += "Parkour queries tried:\n"
    for fq in flowcellqueries:
        msg += f"{fq}\n"

    send_email("Error with flowcell:", msg, config)
    sys.exit("parkour failure.")


def query_parkour_project(config):
    '''
    Given the populated config (including ['data']['projects'])
    give back a dictionary with unique FIDs per project.

    Note that a successful query gives back an object as such:
    {'flowpaths': {'sampleID': ['FID', 'FID2'], 'sampleID2': ['FID3', 'FID4']}'}
    returns a set of flow cell IDs.
    '''
    pid_to_fids = {}
    for project in config['data']['projects']:
        pid = project.split('_')[0]
        res = requests.get(
            config["parkour"]["url"] + f"/api/requests/{pid}/get_flowcell/",
            auth=(config["parkour"]["user"], config["parkour"]["password"]),
            verify=config["parkour"]["pem"]
        )
        if res.status_code == 200:
            parsed = {k: json.loads(v) for k, v in res.json()['flowpaths'].items()}
            pid_to_fids[project] = set([fid for samplel in parsed.values() for fid in samplel])
        else:
            print(f"[red]Parkour query failed for project {project} with status code {res.status_code}[/red]")
            pid_to_fids[project] = set()
    return pid_to_fids