import glob
import os
import shutil
import subprocess as sp

import yaml
from rich import print


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
    perc_used = round(shutil.disk_usage(path).used / shutil.disk_usage(path).total, 2)
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
    

def scan_multiqc(config):
    """
    Collect QC metrices from multiqc json file
    Warning: Assumes 1 project per flowcell
    """
    print("[green]Parse multiqc[/green]")
    QC = {}
    # get json file from multiqc
    qc_dir_root = "Project_" + config["data"]["projects"][0]
    multiqc_path = os.path.join("QC", "multiqc_data", "multiqc_data.json")
    json1 = os.path.join(
        config["info_dict"]["transfer_path"], qc_dir_root, multiqc_path
    )

    if os.path.exists(json1):
        json = json1
    else:
        # for compatibility: revert to old location of QC files
        multiqc_path = os.path.join("multiqc", "multiqc_data", "multiqc_data.json")
        json = os.path.join(
            config["info_dict"]["transfer_path"],
            "FASTQC_" + qc_dir_root,
            multiqc_path,
        )

    if not os.path.exists(json):
        print(f"[red]Warning. json does not exit: {json}[/red]")
        return QC

    jy = yaml.safe_load(open(json))
    dd = jy["report_saved_raw_data"]["multiqc_fastqc"]
    QC["samples"] = list(dd.keys())
    QC["total_sequences"] = [v.get("Total Sequences", None) for v in dd.values()]
    QC["total_bases"] = [v.get("Total Bases", None) for v in dd.values()]

    QC["percent_gc"] = [v.get("%GC", None) for v in dd.values()]
    QC["percent_dedup"] = [
        round(v.get("total_deduplicated_percentage", None), 2) for v in dd.values()
    ]


    screens = list(jy["report_data_sources"].keys())  # list of QC screens (FastQC, ...)
    if "pycoQC" not in screens:
        print(f"[red]Warning. no pycoQC metrics in: {json}[/red]")
    else:
        dd = jy["report_general_stats_data"]["pycoqc"]  # dictionary
        QC["all_median_phred_score"] = [
            round(v.get("all_median_phred_score", None), 2) for v in dd.values()
        ]
        QC["all_median_read_length"] = [
            round(v.get("all_median_read_length", None), 2) for v in dd.values()
        ]
        QC["all_N50"] = [v.get("all_n50", None) for v in dd.values()]

    return QC


def print_header(head):
    left = 50
    right = 50 - len(head)
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