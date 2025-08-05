import glob
import os
import re
import shutil
import subprocess as sp
import sys
from pathlib import Path

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
        'debug': '--debug'
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


def monitor_storage(config):
    """
    collect information on added storage in transfered directory
    and overall storage in partition
    """
    path = config["info_dict"]["transfer_path"]
    SM = {}  # storage monitor

    # get du information for newly added directory
    du_out = sp.check_output(
        ["du", "-sh", path], universal_newlines=True, stderr=sp.STDOUT
    )
    SM["added_storage"] = du_out.split()[0]

    # get partition_info from df: (header, info, '^\n')
    df_lines = (
        sp.check_output(["df", "-BG", path], universal_newlines=True, stderr=sp.STDOUT)
        .strip()
        .split("\n")
    )

    # remove header (= use last line) and split into columns
    df_out = df_lines[-1].split()
    SM["available_storage_GB"], SM["used_storage_perc"] = [df_out[3], df_out[4]]
    return SM


def scan_multiqc(config):
    """
    Collect QC metrices from multiqc json file
    Warning: Assumes 1 project per flowcell
    """
    QC = {}
    # get json file from multiqc
    qc_dir_root = "Project_" + config["data"]["projects"][0]
    multiqc_path = os.path.join("QC", "multiqc", "multiqc_data", "multiqc_data.json")
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

    # get QC metrics from FastQC (pass reads)
    dd = jy["report_saved_raw_data"]["multiqc_fastqc"]
    QC["samples"] = list(dd.keys())
    QC["total_sequences"] = [v.get("Total Sequences", None) for v in dd.values()]
    QC["total_bases"] = [v.get("Total Bases", None) for v in dd.values()]
    # FastQC is the wrong place to get the lengths --> pycoqc
    # QC['sequence_length_range'] = [v.get('Sequence length', None) for v in dd.values()]
    # QC['sequence_length_median'] = [v.get('median_sequence_length', None) for v in dd.values()]
    QC["percent_gc"] = [v.get("%GC", None) for v in dd.values()]
    QC["percent_dedup"] = [
        round(v.get("total_deduplicated_percentage", None), 2) for v in dd.values()
    ]

    # get QC metrics from pycoQC (for all reads)
    # unfortunately jy['report_general_stats_data'] is a list and not a dict
    # first need to get the index of pycoQC (should be 2)
    screens = list(jy["report_data_sources"].keys())  # list of QC screens (FastQC, ...)

    if "pycoQC" not in screens:
        print(f"[red]Warning. no pycoQC metrics in: {json}[/red]")
    else:
        pyco_idx = screens.index("pycoQC")  # index
        dd = jy["report_general_stats_data"][pyco_idx]  # dictionary
        QC["all_median_phred_score"] = [
            round(v.get("all_median_phred_score", None), 2) for v in dd.values()
        ]
        QC["all_median_read_length"] = [
            round(v.get("all_median_read_length", None), 2) for v in dd.values()
        ]
        QC["all_N50"] = [v.get("all_n50", None) for v in dd.values()]
    #       QC['all_bases_pyco'] = [v.get('all_bases', None) for v in dd.values()]
    #       QC['all_reads_pyco'] = [v.get('all_reads', None) for v in dd.values()]

    if "Kraken" not in screens:
        print(f"[red]Warning. no Kraken metrics in: {json}[/red]")
    else:
        kraken_idx = screens.index("Kraken")  # index
        dd = jy["report_general_stats_data"][kraken_idx]  # nested dictionary
        QC["top_species"] = []
        QC["top_percent"] = []

        # loop over all samples (keys) and their dictionaries (values)
        for sample_dict in dd.values():
            sample_items = list(sample_dict.items())
            if len(sample_items) > 0:
                # last element denotes top species and associated percentage
                spec, perc = sample_items[-1]
                # QC['top_species'].append(spec) ## multiQC 'report_general_stats' does not contain species names
                QC["top_percent"].append(round(perc, 2))

        # get top species name (multiQC reports only the top one)
        dd = jy["report_general_stats_headers"][kraken_idx]
        QC["top_species"].append(dd["pct_top_one"]["title"].replace("% ", ""))

    return QC


def fast5_to_pod5(basepath, baseout, cmdlinef):
    """
    searches for 'fast5*' directories in basepath
    runs pod5 conversion
    outputs pod5 combined into baseout/pod5*.
    """
    odir = os.path.join(baseout, "pod5")
    podracercmd = [
        "pod5-convert-from-fast5",
        basepath,  # fast5
        odir,  # pod5out
        "--recursive",
        "--force-overwrite",
        "-p",
        "10",
    ]
    if not os.path.exists("log"):
        os.mkdir("log")
    with open(cmdlinef, "w") as f:
        f.write("#pod5-conversion cmd:\n")
        f.write(" ".join(podracercmd) + "\n")
    sp.check_output(podracercmd)


def merge_pod5(basepath, baseout, cmdlinef):
    """
    searches for 'pod5*' directories in basepath
    runs pod5 merge
    outputs pod5 combined into baseout/pod5*.
    """
    odir = os.path.join(baseout, "pod5")

    # merge pod5 files in both pod5_pass/ and pod5_fail/ (the latter tends to be empty)
    pattern = os.path.join(basepath, "pod5_*", "**", "*.pod5")
    pod5_files = glob.glob(pattern, recursive=True)

    lst_cmd = ["pod5", "merge", pod5_files, os.path.join(odir, "merged.pod5")]

    podracercmd = flatten_irreg_lists(lst_cmd)

    if not os.path.exists("log"):
        os.mkdir("log")
    with open(cmdlinef, "w") as f:
        f.write("#pod5-merge cmd:\n")
        f.write(" ".join(podracercmd) + "\n")
    sp.check_output(podracercmd)


def flatten_irreg_lists(nested_list):
    if isinstance(nested_list, list):
        return [a for i in nested_list for a in flatten_irreg_lists(i)]
    else:
        return [nested_list]


def run_command(cmd, logf):
    """
    run command and print stderr to log file logf
    """
    print(f"[yellow] {cmd} [/yellow]")
    with open(logf, "a") as f:
        f.write("### command: ###\n")
        f.write(cmd + "\n")

    # run
    try:
        with open(logf, "a") as f:
            res = sp.run(cmd, stderr=f, check=True, shell=True)
    except sp.CalledProcessError as e:
        print("[red] Failed with return code [/red]", e.returncode)
        print(res)
        print(f"[red] Check also file {logf} [/red]")
        sys.exit(1)


def guppy2dorado(model_name):
    """
    This is an effort to translate model names
    Currently the model is encoded in the report*.json or can be obtained from
    >pod5 inspect debug <pod5>
    However, those models follow guppy naming conventions and have to be translated to dorado
    Below is a brute force translation which should work for most runs until it fails
    It will ikely have to be extended
    Notice that the there is also a limited set of dorado models that will need update
    """
    new_name = model_name
    new_name = "dna_r10.4.1_e8.2_400bps_sup@v4.1.0"  # most common
    if re.match(r"^rna", model_name):
        new_name = "rna002_70bps_hac@v3"
    elif re.match(r"^dna_r9.4.1", model_name):
        # not the latest model but comaptible with modification calling
        new_name = "dna_r9.4.1_e8_sup@v3.3"
    elif re.match(r"^dna_r10.4.1_e8.2_400bps_5khz_hac_prom", model_name):
        new_name = "dna_r10.4.1_e8.2_400bps_sup@v4.2.0"

    print(f"[red] guppy2dorado: {model_name} -> {new_name}[/red]")
    return new_name


def dorado_basecalling(config, cmdlinef, logf):
    """
    1. get model name
    2. check if basecalling had been started (existing bam file)
    3. run dorado basecalling
    4. create sequencing_summary.txt and put it to fastq/
    5. convert *.bam to *.fastq.gz and put it into fastq/pass
    Notice:
    - dorado does not yet do barcode splitting and only a single bam file will be produced
    - there is no split into pass/ fail/ based on min_qscore (--min_score=0)
    - fastq.gz creation seems an overload, but subsequent workflows depent on it
    """
    # model name = default model (inferred from json)
    model_name = config["info_dict"]["model_def"]

    # try to translate model name to dorado schema
    model_name = guppy2dorado(model_name)

    # overwrite model name if specified in config explicitly
    if config["dorado_basecaller"]["dorado_model"]:
        model_name = config["dorado_basecaller"]["dorado_model"]

    # model directory
    model = os.path.join(config["dorado_basecaller"]["model_directory"], model_name)

    if not os.path.exists(model):
        print(f"[red] Model {model} not available [/red]")
        sys.exit(1)

    # input directory: data
    pod5dir = os.path.join(config["info_dict"]["flowcell_path"], "pod5")

    # output file: currently Dorado will return a single BAM file
    oform = "bam"
    outdir = pod5dir.replace("pod5", oform)
    if not os.path.exists(outdir):
        try:
            os.mkdir(outdir)
        except Exception as e:
            print(f"[red] Creating {outdir}. Error: {e}[/red]")
            sys.exit(1)

    outfile = os.path.join(outdir, "dorado_basecalled." + oform)

    # if the output file (BAM) already exists resume basecalling
    # Notice: resume only works with BAM
    dorado_resume = []
    if oform == "bam" and os.path.exists(outfile):
        print(
            "[yellow] File {} is already available. Resume basecalling [/yellow]".format(
                outfile
            )
        )
        oldfile = outfile.replace(".bam", ".previous.bam")
        shutil.move(outfile, oldfile)
        dorado_resume = ["--resume-from", oldfile]

    # include dorado options passed to config
    dorado_opt = []
    if config["dorado_basecaller"]["dorado_options"]:
        dorado_opt = config["dorado_basecaller"]["dorado_options"].split(" ")

    # discarded effort to have alternative fastq output with '--emit-fastq'
    # if config["dorado_basecaller"]["dorado_output"]=='fastq':
    #    dorado_opt = dorado_opt + ['--emit-fastq']

    cmd = (
        [config["dorado_basecaller"]["dorado_cmd"]]
        + ["basecaller"]
        + dorado_opt
        + dorado_resume
        + [model]
        + [pod5dir]
    )
    cmd.append(f"> {outfile}")
    cmd = " ".join(cmd)
    run_command(cmd, logf)

    # run sequence summary
    # put output into fastq directory (historical reasons: that's where pycoQC will be looking)
    fastq_dir = pod5dir.replace("pod5", "fastq")
    os.makedirs(fastq_dir, exist_ok=True)
    seq_sum = os.path.join(fastq_dir, "sequencing_summary.txt")

    cmd = [config["dorado_basecaller"]["dorado_cmd"], "summary", outfile]
    cmd.append(f"> {seq_sum}")
    cmd = " ".join(cmd)
    run_command(cmd, logf)

    # convert bam to fastq.gz
    if config["dorado_basecaller"]["dorado_output"] == "fastq":
        fastq_dir = os.path.join(fastq_dir, "pass")
        os.makedirs(fastq_dir, exist_ok=True)  # create .../fastq/pass
        bn = os.path.basename(outfile).replace("bam", "fastq.gz")
        fastq_file = os.path.join(fastq_dir, bn)  # .../fastq/pass/*.fastq.gz

        cmd = ["samtools", "bam2fq", outfile, "| gzip >", fastq_file]
        cmd = " ".join(cmd)
        run_command(cmd, logf)


def retRule(rulestr, config):
    return os.path.join(config["paths"]["rulesPath"], rulestr)


def glob2reports(globStr, base_path, flowcell_path):
    reports_dir = os.path.join(flowcell_path, "reports")
    if not os.path.exists(reports_dir):
        print(f"Creating {reports_dir}")
        os.mkdir(reports_dir)
    globber = glob.glob(os.path.join(base_path, globStr))
    # if not globber:
    #     sys.exit("Not {} files found..".format(globStr))
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


def genome_index(config, path):
    """
    Generating genome indices for minimap2
    """
    reference = config["organism"]
    reference_fa = config["genome"][reference]
    cmd = "ln -s " + reference_fa + " " + path + "/" + reference + "_genome.fa;"
    cmd += (
        config["mapping"]["mapping_cmd"]
        + " "
        + config["mapping"]["index_options"]
        + " "
    )
    cmd += path + "/" + reference + "_genome.mmi "
    cmd += path + "/" + reference + "_genome.fa "
    print(cmd)
    sp.check_call(cmd, shell=True)


def overwrite_dir(dir, destbase):
    ret = "copied"
    if os.path.exists(os.path.join(destbase, dir)):
        shutil.rmtree(os.path.join(destbase, dir))
        ret = "replaced"
    shutil.copytree(dir, os.path.join(destbase, dir))
    return ret


def grab_seqsummary(dir):
    globber = glob.glob(os.path.join(dir, "sequencing_summary*"))
    if len(globber) != 1:
        sys.exit("more then 1 sequencing summary.")
    else:
        return globber[0]


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


def config_to_splitseqsummary(config):
    cmd = config["pycoQc"]["barcodeSplit"]
    cmd += " fastq/sequencing_summary.txt "
    cmd += "-o ./"
    return cmd


# deprecated
def config_to_pycoqc(config, seqsum, fqc_sampledir, sampleid, bc_kit):
    cmd = config["pycoQc"]["pycoQc"]
    cmd += " " + seqsum + " "
    cmd += "-o " + os.path.join(fqc_sampledir, "pycoqc_" + sampleid + ".html")
    return cmd


def config_to_mapcmd(config):
    protocol = config["info_dict"]["protocol"]
    org = config["info_dict"]["organism"]

    if protocol in ["dna", "cdna", "rna"] and org in ["mouse", "drosophila", "human"]:
        pref = [config["mapping"]["mapping_cmd"]]
        post = [config["mapping"]["samtools_cmd"]]
        if protocol == "dna":
            pref = pref + config["mapping"]["mapping_dna_options"].split(" ")
            pref.append(config["genome"][config["info_dict"]["organism"]])
            post.append("sort")
            post.append(config["mapping"]["samtools_options"])
            post.append("-o")
            return (pref, post)
        elif protocol == "rna":
            print("protocol = rna")
            pref = pref + config["mapping"]["mapping_rna_options"].split(" ")
            pref = pref + ["--junc-bed", config["transcripts"][org]]
            pref.append(config["genome"][org])
            post.append("sort")

            post.append(config["mapping"]["samtools_options"])
            post.append("-o")
            return (pref, post)
        elif protocol == "cdna":
            print("protocol = cdna")
            pref = pref + config["mapping"]["mapping_rna_options"].split(" ")
            pref = pref + ["--junc-bed", config["transcripts"][org]]
            pref.append(config["genome"][org])
            post.append("sort")
            post.append(config["mapping"]["samtools_options"])
            post.append("-o")
            return (pref, post)
    else:
        return (None, None)
