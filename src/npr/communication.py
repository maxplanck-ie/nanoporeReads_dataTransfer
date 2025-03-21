import glob
import os
import shutil
import smtplib
import sys
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from importlib.metadata import version
from time import sleep

import paramiko
import requests
from rich import print
from scp import SCPClient


def ship_qcreports(config, flowcell):
    """
    Legacy shipment of QC reports. Reconsider since this generates quite some duplication
    The idea is to solve
     1. NGS team can only reads via SAMBA share
     2. bioinfo may not be able to open html on analysis server
    """

    # make shipping dependent on whether sambahost is defined - simplify testing
    # also disables shipping to bioinfocore
    if config["sambahost"]["host"] is None:
        return

    # login info.
    _pkey = paramiko.RSAKey.from_private_key_file(config["sambahost"]["pkey"])
    _user = config["sambahost"]["user"]
    _host = config["sambahost"]["host"]

    # set client.
    client = paramiko.SSHClient()
    policy = paramiko.AutoAddPolicy()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    client.connect(_host, username=_user, pkey=_pkey)
    # Create flowcell folder if it doesn't exist.
    yrstr = flowcell[:4]
    samba_fdir = os.path.join(
        config["paths"]["deepseq_qc"],
        f"Sequence_Quality_{yrstr}",
        f"ONT_{yrstr}",
    )
    # I think latency causes scp to fail so include sleep (legacy from Leily?)
    _cmd = f"mkdir -p {samba_fdir}"
    _stdin, _stdout, _stderr = client.exec_command(_cmd)

    #    print('stdin:')
    #    print(_stdin)
    #    print('stdout:')
    #    print(_stdout)
    print("stderr:")
    print(_stderr)
    sleep(30)

    # copy run_reports & pycoQC
    scp = SCPClient(client.get_transport())
    fc_dir = config["info_dict"]["flowcell_path"]
    fc_name = os.path.basename(fc_dir)
    samba_target = os.path.join(samba_fdir, fc_name)
    bioinfo_target = os.path.join(config["paths"]["bioinfocoredir"], fc_name)

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
        try:
            scp.put(report, samba_target)
        except (paramiko.SSHException, OSError) as e:
            print(f"Error during transfer to sambahost: {str(e)}")

    # copy QC directories: transfer/Project*/QC
    for local_dir in glob.glob(f"{fc_dir}/transfer/Project*/QC/"):
        p_dir = os.path.dirname(os.path.dirname(local_dir))
        p_name = os.path.basename(p_dir)
        target_dir = os.path.join(fc_name, p_name, "QC")

        # copy to bioinfo
        bioinfo_target = os.path.join(config["paths"]["bioinfocoredir"], target_dir)
        print(f"copy {local_dir}")
        shutil.copytree(local_dir, bioinfo_target, dirs_exist_ok=True)

        # copy to sambahost
        samba_target = os.path.join(samba_fdir, target_dir)
        stdin, stdout, stderr = client.exec_command(f"mkdir -p {samba_target}")
        sleep(30)
        if stderr.read():
            print(f"Error creating remote directory {samba_target}. Error: {stderr}")
        try:
            scp.put(local_dir, recursive=True, remote_path=samba_target)
        except (paramiko.SSHException, OSError) as e:
            print(f"Error during transfer of {local_dir} to sambahost: {str(e)}")

    scp.close()
    client.close()


def standard_text(config):
    """
    Draft a short letter to end user that will be part of the success email
    Include also QC metrics obtained from multiqc report
    """
    QC, SM = config["QC"], config["SM"]
    samples = QC.pop("samples")
    frame = "\n-------\n"
    msg = (
        "Dear <...>\n"
        + "The sequencing and first analysis for Project {} is finished\n".format(
            config["data"]["projects"]
        )
        + "\n"
        + "Ouput Folder: {}".format(config["info_dict"]["transfer_path"])
        + frame
        + "Quality Metrics (from mulitQC) \n"
        + f"Samples: {samples}\n"
        + "\n".join([f"{key}: {value}" for key, value in QC.items()])
        + frame
        + "Storage Footprint \n"
        + "\n".join([f"{key}: {value}" for key, value in SM.items()])
        + frame
        + "Please keep in mind that the original pod5 files, needed e.g. to call modified bases, will only be available for the next 6 months and will be permanently deleted afterwards.\n"
        + "Please let me know if something is unclear or if you have any further questions.\n\nKind regards.\n"
    )
    return msg


def send_email(subject, body, config, allreceivers=True):
    """
    Send email including key information about the run
    Also print message to stdout
    """
    mailer = MIMEMultipart("alternative")
    mailer["Subject"] = "[npr] [{}] {} {}".format(
        version("npr"), subject, os.path.basename(config["info_dict"]["base_path"])
    )

    # add standard information from config['data'] and config['info_dict'] to each message
    info = ""
    if "data" in config and "projects" in config["data"]:
        info += "Project: {}\n".format(config["data"]["projects"])
    if "data" in config and "samples" in config["data"]:
        info += "Samples: {}\n".format(config["data"]["samples"])
    if "info_dict" in config:
        info += "\n".join(
            [f"{key}: {value}" for key, value in config["info_dict"].items()]
        )
    if "basecaller" in config:
        info += "\nbasecaller: {}\n".format(config["basecaller"])
    frame = "\n=====\n"
    body = body + frame + info + frame

    mailer["From"] = config["email"]["from"]
    to_email = "to" if allreceivers else "trigger"
    mailer["To"] = config["email"][to_email]
    tomailers = config["email"]["to"].split(",")
    print(f"Email trigger, sending to {tomailers}")
    email = MIMEText(body)
    mailer.attach(email)
    if config["email"]["host"] is not None:
        s = smtplib.SMTP(config["email"]["host"])
        s.sendmail(config["email"]["from"], tomailers, mailer.as_string())
    print("[green]Subject: {}\n{}[/green]".format(mailer["Subject"], body))


def query_parkour(config, flowcell, msg):
    """
    query parkour to update info_dict (protocol and organism)
    """

    if not config["parkour"]["url"]:
        msg += "Parkour URL not specified."
        return msg

    # Old manual escapes.
    if flowcell == "20221014_1045_X5_FAV39027_f348bc5c":
        fc = "FAV39027_reuse"
    elif flowcell == "20221107_1020_X3_FAV08360_71e3fa80":
        fc = "FAV08360-1"
    elif flowcell == "20230331_1220_X4_FAV22714_872a401d":
        fc = "FAV22714-2"
    else:
        try:
            fc = flowcell.split("_")[3]
        except Exception as e:
            print(f"[red]Exception: {e}[/red]")
            msg += 'Parkour Error: flowcell "{}" cannot be queried in parkour.'.format(
                flowcell
            )
            send_email("Error with flowcell:", msg, config, allreceivers=False)
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
            info_dict = {}
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
            organism = parkour_dict[first_key][first_entry][3][1]
            protocol = parkour_dict[first_key][first_entry][1]

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
                print("protocol not found. Default to dna.")
                protocol = "dna"

            # update config['info_dict']
            config["info_dict"]["protocol"] = protocol
            if organism not in config["genome"].keys():
                organism = "other"
            config["info_dict"]["organism"] = str(organism)
            return msg

    msg += f"Parkour query failed for {fc}.\n"
    msg += "Parkour queries tried:\n"
    for fq in flowcellqueries:
        msg += f"{fq}\n"

    send_email("Error with flowcell:", msg, config)
    sys.exit("parkour failure.")
