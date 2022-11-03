import os
import shutil
import sys
from rich import print
import glob

def retRule(rulestr, config):
    return(
        os.path.join(
            config['paths']['rulesPath'],
            rulestr
        )
    )

def glob2reports(globStr, base_path, flowcell_path):
    reports_dir = os.path.join(
        flowcell_path,
        'reports'
    )
    if not os.path.exists(reports_dir):
        print("Creating {}".format(reports_dir))
        os.mkdir(reports_dir)
    globber = glob.glob(
            os.path.join(
                base_path,
                globStr
            )
    )
    if not globber:
        sys.exit("Not {} files found..".format(globStr))
    for s in globber:
        dest = os.path.join(
            flowcell_path,
            'reports',
            os.path.basename(s)
        )
        if not os.path.exists(dest):
            shutil.copy(
                s,
                dest
                )

def get_seqdir(groupdir, seqdir):
    max_num = 0
    for dir in glob.glob(os.path.join(groupdir, seqdir+"*")):
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
    cmd += config["mapping"]["mapping_cmd"] +" "+ config["mapping"]["index_options"] + " "
    cmd += path + "/" + reference + "_genome.mmi "
    cmd += path + "/" + reference + "_genome.fa "
    print(cmd)
    sp.check_call(cmd, shell=True)

def overwrite_dir(dir, destbase):
    ret = 'copied'
    if os.path.exists(
        os.path.join(destbase, dir)
    ):
        shutil.rmtree(
            os.path.join(destbase, dir)
        )
        ret = 'replaced'
    shutil.copytree(
        dir,
        os.path.join(destbase, dir)
    )
    return(ret)

def grab_seqsummary(dir):
    globber = glob.glob(
        os.path.join(
            dir,
            'sequencing_summary*'
        )
    )
    if len(globber) != 1:
        sys.exit("more then 1 sequencing summary.")
    else:
        return(globber[0])

def getfoot(dir):
    footprint = 0
    for path, dirs, files in os.walk(dir):
        for f in files:
            footprint += os.path.getsize(
                os.path.join(
                    path,f
                )
            )
    return(footprint)

# return commands.
def config_to_basecallcmd(config):
    # base cmd.
    cmd = config['guppy_basecaller']['base_calling_cmd']
    # in - and output folder.
    cmd += " -i {}".format(
        config["info_dict"]["poddir"]
    )
    cmd += " -s fastq"
    barcode = False if config['bc_kit'] == 'no_bc' else True
    cmd += " -c {} ".format(config['info_dict']['model'])
    cmd += config["guppy_basecaller"]["base_calling_options"]
    if barcode is True:
        cmd += " --barcode_kits {} ".format(config["bc_kit"])
        cmd += config["guppy_basecaller"]["base_calling_barcode_options"]
    else:
        cmd += " --trim_strategy dna "
    return (cmd)

def config_to_splitseqsummary(config):
    cmd = config["pycoQc"]["barcodeSplit"]
    cmd += " fastq/sequencing_summary.txt "
    cmd += "-o ./"
    return(cmd)

def config_to_pycoqc(config, seqsum, fqc_sampledir, sampleid, bc_kit):
    cmd = config["pycoQc"]["pycoQc"]
    cmd += ' ' + seqsum + ' '
    cmd += "-o " + os.path.join(
        fqc_sampledir,
        'pycoqc_' + sampleid + '.html'
    )
    return(cmd)

def config_to_mapcmd(config):
    protocol = config['info_dict']['protocol']
    org = config['info_dict']['organism']

    if protocol in [
        'dna',
        'cdna',
        'rna'
    ] and org in [
        'mouse',
        'drosophila',
        'human'
    ] :
        pref = [config['mapping']['mapping_cmd']]
        post = [config['mapping']['samtools_cmd']]
        if protocol == 'dna':
            pref = pref + config['mapping']['mapping_dna_options'].split(' ')
            pref.append(config['genome'][config['info_dict']['organism']])
            post.append('sort')
            post.append(
                config['mapping']['samtools_options']
            )
            post.append(
                '-o'
            )
            return(pref,post)
        elif protocol == 'rna':
            print('protocol = rna')
            pref = pref + config['mapping']['mapping_rna_options'].split(' ')
            pref = pref + ['--junc-bed', config['transcripts'][org]]
            pref.append(config['genome'][org])
            post.append('sort')

            post.append(
                config['mapping']['samtools_options']
            )
            post.append(
                '-o'
            )
            return(pref,post)
        elif protocol == 'cdna':
            print('protocol = cdna')
            pref = pref + config['mapping']['mapping_rna_options'].split(' ')
            pref = pref + ['--junc-bed', config['transcripts'][org]]
            pref.append(config['genome'][org])
            post.append('sort')
            post.append(
                config['mapping']['samtools_options']
            )
            post.append(
                '-o'
            )
            return(pref,post)
    else:
        return(None,None)

