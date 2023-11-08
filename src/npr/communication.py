from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib
from rich import print
import requests
import sys
import os
import glob
import shutil
from importlib.metadata import version
import paramiko
from scp import SCPClient
from time import sleep

def ship_qcreports(config, flowcell):
    '''
    update samba drive with a flowcell folder,
    copy pycoQC reports and run report into samba drive
    copy both reports in bioinfo qc drive as well
    '''

    # make shipping dependent on whether sambahost is defined - simplify testing
    # also disables shipping to bioinfocore
    if config['sambahost']['host'] is None:
        return

    # login info.
    _pkey = paramiko.RSAKey.from_private_key_file(
        config['sambahost']['pkey']
    )
    _user = config['sambahost']['user']
    _host = config['sambahost']['host']
    
    # set client.
    client = paramiko.SSHClient()
    policy = paramiko.AutoAddPolicy()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    
    client.connect(
        _host,
        username=_user,
        pkey=_pkey
    )
    # Create flowcell folder if it doesn't exist.
    yrstr = flowcell[:4]
    samba_fdir = os.path.join(
        config['paths']['deepseq_qc'],
        'Sequence_Quality_{}'.format(yrstr),
        'ONT_{}'.format(yrstr),
        flowcell
    )
    # I think latency causes scp to fail so include sleep
    _cmd = 'mkdir -p {}'.format(samba_fdir)
    _stdin, _stdout, _stderr = client.exec_command(_cmd)
    print('stdin:')
    print(_stdin)
    print('stdout:')
    print(_stdout)
    print('stderr:')
    print(_stderr)
    sleep(30)

    # copy run_reports & pycoQC
    scp = SCPClient(client.get_transport())

    pycoqcs = glob.glob(
        os.path.join(config['info_dict']['flowcell_path'], 'FASTQC*', '*', '*pycoqc.html')
    ) + glob.glob(
        os.path.join(config['info_dict']['flowcell_path'], 'reports', '*.html')
    )
    for qcreport in pycoqcs:
        basename = os.path.basename(qcreport)
        if 'pycoqc' in basename:
            sampleID = qcreport.split('/')[-2].replace('Sample_', '')
            projID = 'Project_' + qcreport.split('/')[-3].split('_')[2]
            basename = projID + '_' + sampleID + '_' + basename
        sambadest = os.path.join(samba_fdir, basename)
        print('Trying to copy {} to {}'.format(qcreport, sambadest))
        bioinfodest = os.path.join(config['paths']['bioinfocoredir'], flowcell + '_' + basename)
        scp.put(qcreport, sambadest)
        shutil.copyfile(qcreport, bioinfodest)
    scp.close()
    client.close()


def standard_text(config):
    """
    Draft a short letter to end user that will be part of the success email
    Include also QC metrics obtained from multiqc report
    """
    QC, SM = config['QC'], config['SM']
    samples = QC.pop('samples')
    frame = "\n-------\n"
    msg="Dear <...>\n" +\
     "The sequencing and first analysis for Project {} is finished\n".format(config['data']['projects']) +\
     "\n" +\
    "Ouput Folder: {}".format(config['info_dict']['transfer_path']) +\
    frame +\
    "Quality Metrics (from mulitQC) \n" +\
    "Samples: {}\n".format(samples) +\
    "\n".join([f"{key}: {value}" for key, value in QC.items()]) +\
    frame +\
    "Storage Footprint \n" +\
    "\n".join([f"{key}: {value}" for key, value in SM.items()]) +\
    frame +\
    "Please let me know if something is unclear or if you have any other questions.\n\nKind regards.\n"
    return msg 


def send_email(subject, body, config, allreceivers=True):
    """
    Send email including key information about the run
    Also print message to stdout
    """
    mailer = MIMEMultipart('alternative')
    mailer['Subject'] = "[npr] [{}] {} {}".format(
        version('npr'),
        subject,
        os.path.basename(config['info_dict']['base_path'])
    )

    # add standard information from config['data'] and config['info_dict'] to each message
    info = ""
    if 'data' in config and 'projects' in config['data']:
        info += 'Project: {}\n'.format(config['data']['projects'])
    if 'data' in config and 'samples' in config['data']:
        info += 'Samples: {}\n'.format(config['data']['samples'])
    if 'info_dict' in config:
        info += "\n".join([f"{key}: {value}" for key, value in config['info_dict'].items()])
    if 'basecaller' in config:
        info += '\nbasecaller: {}\n'.format(config['basecaller'])
    frame = "\n=====\n"
    body = body + frame + info + frame

    mailer['From'] = config['email']['from']
    to_email = 'to' if allreceivers else 'trigger'
    mailer['To'] = config['email'][to_email]
    tomailers = config['email']['to'].split(',')
    print("Email trigger, sending to {}".format(tomailers))
    email = MIMEText(body)
    mailer.attach(email)
    if config['email']['host'] is not None:
        s = smtplib.SMTP(config['email']['host'])
        s.sendmail(
            config['email']['from'],
            tomailers,
            mailer.as_string()
        )
    print('[green]Subject: {}\n{}[/green]'.format(mailer['Subject'], body)) 


def query_parkour(config, flowcell, msg):
    """
    query parkour to update info_dict (protocol and organism)
    """
    if not config["parkour"]["url"]:
        msg += "Parkour URL not specified."
        return(msg)

    # Old manual escapes.
    if flowcell == '20221014_1045_X5_FAV39027_f348bc5c':
        fc = 'FAV39027_reuse'
    if flowcell == '20221107_1020_X3_FAV08360_71e3fa80':
        fc = 'FAV08360-1'
    if flowcell == '20230331_1220_X4_FAV22714_872a401d':
        fc = 'FAV22714-2'
    else:
        try:
            fc = flowcell.split("_")[3]
        except Exception as e:
            print('[red]Exception: {}[/red]'.format(e))
            msg += 'Parkour Error: flowcell "{}" cannot be queried in parkour.'.format(flowcell)
            send_email('Error with flowcell:', msg, config)
            sys.exit(1)
    # test for flow cell re-use.
    # flowcell that's re-used gets higher increment.
    postfixes = ['-5', '-4', '-3', '-2', '-1','']
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
            if fc == 'PAK78871' or fc == 'PAK79330' or fc =='PAK77043' or fc == 'FAV22714-2' or fc == 'PAK78965' or fc == 'PAK79757':
                print("[red] Protocol override [/red]")
                protocol = 'cdna'
            if 'cDNA' in protocol:
                protocol = 'cdna'
            elif 'DNA' in protocol:
                protocol = 'dna'
            elif "RNA" in protocol:
                protocol = 'rna'
            else:
                print('protocol not found. Default to dna.')
                protocol = 'dna'

            # update config['info_dict']
            config['info_dict']['protocol'] = protocol
            if organism not in config['genome'].keys():
                organism = "other"
            config['info_dict']['organism'] = str(organism)
            return (msg)

    msg += "Parkour query failed for {}.\n".format(fc)
    msg += "Parkour queries tried:\n"
    for fq in flowcellqueries:
        msg += "{}\n".format(fq)
    send_email('Error with flowcell:', msg, config)
    sys.exit("parkour failure.")
