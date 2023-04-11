from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib
from rich import print
import requests
import sys
import os
from importlib.metadata import version


def send_email(body, version, flowcell, config, allreceivers=True):
    mailer = MIMEMultipart('alternative')
    mailer['Subject'] = "[npr] [{}] {}".format(
        version,
        flowcell
    )
    mailer['From'] = config['email']['from']
    to_email = 'to' if allreceivers else 'trigger'
    mailer['To'] = config['email'][to_email]
    tomailers = config['email']['to'].split(',')
    print("Email trigger, sending to {}".format(tomailers))
    email = MIMEText(body)
    mailer.attach(email)
    s = smtplib.SMTP(config['email']['host'])
    s.sendmail(
        config['email']['from'],
        tomailers,
        mailer.as_string()
    )

def query_parkour(config, flowcell, msg):
    """
    query parkour.
    """
    # Old manual escapes.
    if flowcell == '20221014_1045_X5_FAV39027_f348bc5c':
        fc = 'FAV39027_reuse'
    if flowcell == '20221107_1020_X3_FAV08360_71e3fa80':
        fc = 'FAV08360-1'
    if flowcell == '20230331_1220_X4_FAV22714_872a401d':
        fc = 'FAV22714-2'
    else:
        fc = flowcell.split("_")[3]
    # test for flow cell re-use.
    # flowcell that's re-used gets higher increment.
    postfixes = ['-5', '-4', '-3', '-2', '-1']
    flowcellqueries = []
    for pf in postfixes:
        d = {'flowcell_id': fc + pf}
        flowcellqueries.append(fc+pf)
        res = requests.get(
            config["parkour"]["url"],
            auth=(
                config["parkour"]["user"],
                config["parkour"]["password"]
            ),
            params=d,
            verify=config['parkour']['pem']
        )
        if res.status_code == 200:
            info_dict = {}
            msg += "Parkour query 200 for fid {} with re-use index {}.\n".format(fc, pf)
            msg += "\n"
            msg += "Parkour queried with:\n"
            for fq in flowcellqueries:
                msg += "{}\n".format(fq)
            msg += "\n"
            parkour_dict = res.json()
            print(parkour_dict)
            first_key = list(parkour_dict.keys())[0]
            first_entry = list(parkour_dict[first_key].keys())[0]
            organism = parkour_dict[first_key][first_entry][-3]
            protocol = parkour_dict[first_key][first_entry][1]
            if fc == 'PAK78871' or fc == 'PAK79330' or fc =='PAK77043' or fc == 'FAV22714-2':
                print("[red] Protocol override [/red]")
                protocol = 'cdna'
            if 'cDNA' in protocol:
                protocol = 'cdna'
            elif 'DNA' in protocol:
                protocol = 'dna'
            elif "RNA" in protocol:
                protocol = 'rna'
            else:
                print('protocol not found Default to dna.')
                protocol = 'dna'
                #sys.exit("protocol not found")
            info_dict['protocol'] = protocol
            if organism not in config['genome'].keys():
                organism = "other"
            info_dict["organism"] = str(organism)
            info_dict["protocol"] = protocol
            msg += "nucleic acid  type = {}\n".format(protocol)
            msg += "organism = {}\n".format(organism)
            return (info_dict, msg, False)
    info_dict = {}
    msg += "Parkour query failed for {}.\n".format(fc)
    msg += "Parkour queries tried:\n"
    for fq in flowcellqueries:
        msg += "{}\n".format(fq)
    send_email(
        msg,
        version('npr'),
        flowcell,
        config
    )
    return (info_dict, msg, True)
    sys.exit("parkour failure.")