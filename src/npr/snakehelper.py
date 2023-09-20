import os
import shutil
import sys
from rich import print
import glob
import subprocess as sp
from pathlib import Path

def fast5_to_pod5(basepath, baseout, cmdlinef):
    '''
    searches for 'fast5*' directories in basepath
    runs pod5 conversion 
    outputs pod5 combined into baseout/pod5*.
    '''
    odir = os.path.join(
        baseout,
        'pod5'
    )
    podracercmd = [
        'pod5-convert-from-fast5',
        basepath, # fast5
        odir, # pod5out
        '--recursive',
        '--force-overwrite',
        '-p',
        '10'
    ]
    if not os.path.exists('log'):
        os.mkdir('log')
    with open(cmdlinef, 'w') as f:
        f.write('#pod5-conversion cmd:\n')
        f.write(' '.join(podracercmd) + '\n')
    sp.check_output(podracercmd)



def merge_pod5(basepath, baseout, cmdlinef):
    '''
    searches for 'pod5*' directories in basepath
    runs pod5 merge 
    outputs pod5 combined into baseout/pod5*.
    '''
    odir = os.path.join(baseout, 'pod5')

    pod5_files=glob.glob(os.path.join(basepath, "pod5_pass", "*.pod5"))
    
    lst_cmd = [
        'pod5', 'merge',
        pod5_files,
        os.path.join(odir, "merged.pod5")  
    ]
    
    podracercmd = flatten_irreg_lists(lst_cmd)

    if not os.path.exists('log'):
        os.mkdir('log')
    with open(cmdlinef, 'w') as f:
        f.write('#pod5-merge cmd:\n')
        f.write(' '.join(podracercmd) + '\n')
    sp.check_output(podracercmd)


def flatten_irreg_lists(nested_list):
    if isinstance(nested_list, list):
        return [a for i in nested_list for a in flatten_irreg_lists(i)]
    else:
        return [nested_list]


def dorado_basecalling(config, cmdlinef, logf):

    # model directory
    model=os.path.join(
        config["dorado_basecaller"]["model_directory"],
        config["dorado_basecaller"]["dorado_model"]
    )

    if not os.path.exists(model):
        print("[red] Model {} not available [/red]".format(model))
        sys.exit(1)

    # input directory: data
    pod5dir = os.path.join(
        config['info_dict']['flowcell_path'],
        'pod5'
    )

    # output file: bam or fastq
    oform = config["dorado_basecaller"]["dorado_output"]  # 'fastq' or 'bam'
    outdir = pod5dir.replace('pod5', oform)
    if not os.path.exists(outdir):
        try:
            os.mkdir(outdir)
        except Exception as e:
            print("[red] Creating {}. Error: {}[/red]".format(outdir,e))
            sys.exit(1)

    outfile =  os.path.join(
        outdir,
        'dorado_basecalled.' + oform
    )

    # if the output file (BAM) already exists resume basecalling
    # Notice: resume only works with BAM
    dorado_resume = []
    if oform == 'bam' and os.path.exists(outfile):
        print("[yellow] File {} is already available. Resume basecalling [/yellow]".format(outfile))
        oldfile = outfile.replace('.bam' , '.previous.bam')
        shutil.move(outfile, oldfile)
        dorado_resume = [ '--resume-from' , oldfile ]

 
    # include dorado options passed to config
    dorado_opt= []
    if config["dorado_basecaller"]["dorado_options"]:
        dorado_opt = config["dorado_basecaller"]["dorado_options"].split(' ')

    if config["dorado_basecaller"]["dorado_output"]=='fastq':
        dorado_opt = dorado_opt + ['--emit-fastq']

    cmd = [ config['dorado_basecaller']['dorado_cmd'] ] +\
        [ 'basecaller' ] +\
        dorado_opt +\
        dorado_resume +\
        [ model ] +\
        [ pod5dir ]
    breakpoint()
    
    # If output is fastq also compress
    if config["dorado_basecaller"]["dorado_output"]=='fastq':
        cmd.append('| gzip > {}.gz'.format(outfile))
    else:
        cmd.append('> {}'.format(outfile))

    cmd = ' '.join(cmd)
    print('[yellow] {} [/yellow]'.format(cmd))
    with open(logf, 'a') as f:
        f.write("#dorado-basecaller cmd:\n")
        f.write(cmd + '\n')

    try:
        with open(logf, "a") as f:
            res = sp.run(cmd, stderr=f, shell=True)
    except sp.CalledProcessError as e:
        print("[red] Dorado failed with return code [/red]", e.returncode)
        print(res)
        print("[red] Check also file {} [/red]".format(logf))
        sys.exit(1)



def guppy_basecalling(config, cmdlinef, logf):
    pod5dir = os.path.join(
        config['info_dict']['flowcell_path'],
        'pod5'
    )
    cmd = [config['guppy_basecaller']['base_calling_cmd']] +\
        [
        '-i',
        pod5dir,
        '-s',
        pod5dir.replace('pod5', 'fastq'),
        '-c',
        config['info_dict']['model']
        ] +\
        config["guppy_basecaller"]["base_calling_options"].split(' ')
    barcode = False if config['bc_kit'] == 'no_bc' else True
    if barcode is True:
        cmd = cmd +\
            [
            '--barcode_kits',
            config['bc_kit']
            ] +\
            config["guppy_basecaller"]["base_calling_barcode_options"].split(' ')
    else:
        cmd = cmd +\
            [
            '--trim_strategy',
            'dna'
            ]
    if config['info_dict']['protocol'] == 'rna':
        cmd = cmd +\
            config['guppy_basecaller']['base_calling_RNA_options'].split(' ')
    with open(logf, 'a') as f:
        f.write("#guppy-basecaller cmd:\n")
        f.write(' '.join(cmd) + '\n')
    print(cmd)
    sp.check_output(cmd)

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
    # if not globber:
    #     sys.exit("Not {} files found..".format(globStr))
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


def getfast5foot(f5dir, pod5dir):
    if os.path.exists(os.path.join(f5dir, 'fast5*')):
        f5foot = get_size_of_files(f5dir, 'fast5*')
    if os.path.exists(os.path.join(pod5dir, 'pod5*')):
        pod5foot = get_size_of_files(pod5dir, 'pod5*')
        try:
            return(round(pod5foot/f5foot, 2))
        except ZeroDivisionError:
            return(float('nan'))

def get_size_of_files(in_dir, dir_name):
    size_files = 0
    for d in glob.glob(os.path.join(in_dir, dir_name)):
        for path, dirs, files in os.walk(d):
            for f in files:
                size_files += os.path.getsize(
                    os.path.join(path,f)
                )
    return(size_files)


def config_to_splitseqsummary(config):
    cmd = config["pycoQc"]["barcodeSplit"]
    cmd += " fastq/sequencing_summary.txt "
    cmd += "-o ./"
    return(cmd)

# deprecated
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

