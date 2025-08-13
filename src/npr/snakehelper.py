import glob
import os
import shutil
from pathlib import Path
import pandas as pd
import json
from rich import print
import datetime

def config_to_smkcmd(snakemake_config):
    args = []
    
    # Bool flags and their corresponding CLI options.
    bool_flags = {
        'dryrun': '--dryrun',
        'verbose': '--verbose', 
        'printshellcmds': '--printshellcmds',
        'use_conda': '--use-conda',
        'printdag': '--printdag',
        'debug': '--debug',
        'rerun-incomplete': '--rerun-incomplete'
    }
    
    value_options = {
        'cores': '-c',
        'conda_prefix': '--conda-prefix',
        'max_jobs_per_second': '--max-jobs-per-second'
    }
    
    for key, value in snakemake_config.items():
        if key == 'snakefile':
            args.extend(['-s', value])
        elif key in bool_flags and value:
            args.append(bool_flags[key])
        elif key in value_options:
            args.extend([value_options[key], str(value)])
        elif key == 'rerun_triggers' and isinstance(value, list):
            for trigger in value:
                args.extend(['--rerun-triggers', trigger])
    
    return args


def get_disk_stat(path, checkprint = False):
    perc_used = round(shutil.disk_usage(path).used / shutil.disk_usage(path).total, 2) * 100
    free = round(shutil.disk_usage(path).free / (1024 ** 3), 2)
    used = round(shutil.disk_usage(path).used / (1024 ** 3), 2)
    if checkprint:
        added = round(sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(path) for f in fs) / (1024**3), 2)
        return(
            {
                "percentage": perc_used,
                "free_GB": free,
                "used_GB": used,
                "added_GB": added,
            }
        )
    return(
        {
            "percentage": perc_used,
            "free_GB": free,
            "used_GB": used,
        }
    )
    

def get_qc(config):
    """
    Collect QC metrices from a flow cell.
    Idea is to get (per sample)
     - total number of reads
     - N50
     - total bps
     - median read length
     - median Q score
     - kraken organism vs parkour organism
    
    And include more general information:
     - model
     - parkour protocol
     - modifications
     - reference genome
     - flowcell
     - kit
     - barcoding
     - barcode_kit
    """
    print("[green]get_qc[/green]")

    QC = {}
    # GI will contain the general information. Single layer information
    # QC will contain sample specific metrics, so lists per parameter.
    QC['GI'] = {}
    # Metrics
    # General information
    QC['GI']['Model'] = config['info_dict']['model']
    QC['GI']['Parkour Protocol'] = config['info_dict']['parkour_protocol']
    QC['GI']['Modifications'] = config['info_dict']['modifications']
    QC['GI']['Reference genome'] = config['info_dict']['organism_label']
    QC['GI']['Flowcell'] = config['info_dict']['flowcell']
    QC['GI']['Kit'] = config['info_dict']['kit']
    QC['GI']['Barcoding'] = config['info_dict']['barcoding']
    QC['GI']['Barcode Kit'] = config['info_dict']['barcode_kit']

    # Initiate GC for all samples
    QC['QC'] = {}
    for sid in config['data']['samples']:
        QC['QC'][sid] = {
            'Sample name': config['data'][sid]['Sample_Name'],
            'Project': config['data'][sid]['Sample_Project'],
            'Kraken top hit': None,
            'Parkour organism': config['info_dict']['organism'],
            'Total bp': None,
            'Total reads': None,
            'N50': None,
            'Median length': None,
            'Median Q': None,
            '% reads Q >= 18': '0%'
        }
    for project in config['data']['projects']:
        QCpath = Path(config['info_dict']['transfer_path'], f'Project_{project}', 'QC')
        # kraken
        for kraken_report in QCpath.glob("*kraken.report"):
            _sampleid = kraken_report.stem.split('_')[0]
            assert _sampleid in QC['QC'], f"Sample ID {_sampleid} not found in QC['QC']"
            krakdf = pd.read_csv(kraken_report, sep='\t', header=None)
            QC['QC'][_sampleid]['Kraken top hit'] = krakdf.iloc[krakdf[2].idxmax()][5].replace(' ', '')
        # pycoqc metrics
        for pycoqc_report in QCpath.glob("*pycoqc.json"):
            _sampleid = pycoqc_report.stem.split('_')[0]
            assert _sampleid in QC['QC'], f"Sample ID {_sampleid} not found in QC['QC']"
            with open(pycoqc_report) as f:
                pycoqc_data = json.load(f)
                QC['QC'][_sampleid]['Total reads'] = pycoqc_data['All Reads']['basecall']['reads_number']
                QC['QC'][_sampleid]['Total bp'] = pycoqc_data['All Reads']['basecall']['bases_number']
                QC['QC'][_sampleid]['N50'] = pycoqc_data['All Reads']['basecall']['N50']
                QC['QC'][_sampleid]['Median length'] = pycoqc_data['All Reads']['basecall']['len_percentiles'][49]
                QC['QC'][_sampleid]['Median Q'] = pycoqc_data['All Reads']['basecall']['qual_score_percentiles'][49]
                for index, qscore in enumerate(pycoqc_data['All Reads']['basecall']['qual_score_percentiles']):
                    if qscore >= 18:
                        QC['QC'][_sampleid]['% reads Q >= 18'] = f'{index}%'
                        break
    _elapsed = datetime.datetime.now() - config['start_time']
    days = _elapsed.days
    hours, remainder = divmod(_elapsed.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    QC['GI']['Elapsed time'] = f"{days}d {hours}h {minutes}m {seconds}s"
    return QC


def print_header(head):
    left = 40
    right = 40 - len(head)
    print(f"\n\n[bold] {left * '-'} {head} {right * '-'} [/bold]")

def glob2reports(globStr, base_path, flowcell_path):
    reports_dir = os.path.join(flowcell_path, "reports")
    if not os.path.exists(reports_dir):
        print(f"Creating {reports_dir}")
        os.mkdir(reports_dir)
    globber = glob.glob(os.path.join(base_path, globStr))
    for s in globber:
        dest = os.path.join(flowcell_path, "reports", os.path.basename(s))
        if not os.path.exists(dest):
            shutil.copy(s, dest)


def get_seqdir(groupdir, seqdir):
    """
    Function returns proper target directory if a group has multiple
    seqdir folders (sequencing_data1,2,3...)
    """
    max_num = 0
    for dir in glob.glob(os.path.join(groupdir, seqdir + "*")):
        if dir.split(seqdir)[-1] != "":
            num = dir.split(seqdir)[-1]
            max_num = max(max_num, int(num))
    if max_num != 0:
        seqdir = seqdir + str(max_num)
    if not os.path.exists(os.path.join(groupdir, seqdir, "OxfordNanopore")):
        os.makedirs(os.path.join(groupdir, seqdir, "OxfordNanopore"))
    return os.path.join(groupdir, seqdir, "OxfordNanopore")


def getfast5foot(f5dir, pod5dir):
    if os.path.exists(os.path.join(f5dir, "fast5*")):
        f5foot = get_size_of_files(f5dir, "fast5*")
    if os.path.exists(os.path.join(pod5dir, "pod5*")):
        pod5foot = get_size_of_files(pod5dir, "pod5*")
        try:
            return round(pod5foot / f5foot, 2)
        except ZeroDivisionError:
            return float("nan")


def get_size_of_files(in_dir, dir_name):
    size_files = 0
    for d in glob.glob(os.path.join(in_dir, dir_name)):
        for path, dirs, files in os.walk(d):
            for f in files:
                size_files += os.path.getsize(os.path.join(path, f))
    return size_files