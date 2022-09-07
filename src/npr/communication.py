from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib
from rich import print
import requests
import sys
import os
from importlib.metadata import version


def send_email(body, version, flowcell, config):
    mailer = MIMEMultipart('alternative')
    mailer['Subject'] = "[npr] [{}] {}".format(
        version,
        flowcell
    )
    mailer['From'] = config['email']['from']
    mailer['To'] = config['email']['to']
    email = MIMEText(body)
    mailer.attach(email)
    s = smtplib.SMTP(config['email']['host'])
    s.sendmail(
        config['email']['from'],
        config['email']['to'].split(','),
        mailer.as_string()
    )

def query_parkour(config, flowcell, msg):
    """
    query parkour.
    """
    fc = flowcell.split("_")[3]
    d = {'flowcell_id': fc}
    res = requests.get(
        config["parkour"]["url"],
        auth=(
            config["parkour"]["user"],
            config["parkour"]["password"]
        ),
        params=d,
        verify=config['parkour']['pem']
    )

    if res.status_code == 200: # Flowcell exists!
        info_dict = {}
        msg += "Parkour query 200.\n"
        parkour_dict = res.json()
        first_key = list(parkour_dict.keys())[0]
        first_entry = list(parkour_dict[first_key].keys())[0]
        organism = parkour_dict[first_key][first_entry][-3]
        protocol = parkour_dict[first_key][first_entry][1]
        if 'cDNA' in protocol:
            protocol = 'cdna'
        elif 'DNA' in protocol:
            protocol = 'dna'
        elif "RNA" in protocol:
            protocol = 'rna'
        else:
            sys.exit("protocol not found")
        info_dict['protocol'] = protocol
        if organism not in config['genome'].keys():
            organism = "other"
        info_dict["organism"] = str(organism)
        info_dict["protocol"] = protocol
        return (info_dict, msg)
    else:
        msg += "Parkour query failed. Crash'n'Burn.\n"
        send_email(
            msg,
            version('npr'),
            flowcell,
            config
        )
        sys.exit("parkour failure.")