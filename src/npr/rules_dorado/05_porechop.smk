'''
Adaptor trimming with porechop and summary
'''

# define source and target pattern
source = sample_dat + ".fastq.gz"
subset_fastq = sample_dat + "_subset.fastq"
target_fastq = sample_dat + "_porechop.fastq"
target_info  =  sample_qc + "_porechop.info"
logpat = sample_log + "_porechop.log"
bchpat = sample_bch + "_porechop.tsv"

# Below are some functions to extract more information from the porechop.info file (requires -v 1 (default))
# Notice that Porechop output is rather unstructured
# The function below create files that are compatabile with multiqc (see: multiqc_config.yaml)
def get_trimmed(filename,info):
    '''
    get the list of trimmed barcodes, adaptors and their associated sequence
    Trimming adapters from read ends
        SQK-MAP006_Y_Top_SK63: GGTTGTTTCTGTTGGTGCTGATATTGCT
        SQK-MAP006_Y_Bottom_SK64: GCAATATCAGCACCAACAGAAA
    '''
    suffix = '_porechop.info'
    sample_name = os.path.basename(filename).rsplit(suffix)[0]
    pattern = re.compile(r'^(.*?):(.*?)$')

    info['trimmed'][sample_name]={}
    with open(filename, 'r') as input_file:
        for line in input_file:
            line  = re.sub(r'\x1b\[[0-9;]*m', '', line).strip()
            match = pattern.match(line)
            if match:
                seq_name = match.group(1).strip()
                seq_value = match.group(2).strip()
                info['trimmed'][sample_name][seq_name]=seq_value


def get_percent(filename,info):
    '''
    Read maximal percentage identities for all adaptors and barcodes
    Notice that Porechop formatting contains wild spaces and annyoing ANSI escapes
    A typical line reads:
        "PCR Adaptor 1"    "85.5"   "79.0"
    '''
    suffix = '_porechop.info'
    sample_name = os.path.basename(filename).rsplit(suffix)[0]
    info['best_start'][sample_name]={}
    info['best_end'][sample_name]={}
    pattern = re.compile(r'^(.*?)\s+(\d+\.\d+)\s+(\d+\.\d+)$')

    with open(filename, 'r', newline='') as input_file:
        for line in input_file:
            # Clean ANSI escapes (colours, underlines, ...)
            line  = re.sub(r'\x1b\[[0-9;]*m', '', line).strip()
            match = pattern.match(line)
            if match:
                value1 = match.group(2)
                value2 = match.group(3)
                seq_name = match.group(1).strip()
                # replace internal white spaces by "_"
                seq_name = "_".join(seq_name.split())
                info['best_start'][sample_name][seq_name] = value1
                info['best_end'][sample_name][seq_name] = value2


def print_info(info, dir_name, key='trimmed'):
    '''
    print diagnostic vector for each samples 1 sample = 1 row
    key = 'trimmed' (default), 'best_start' or 'best_end'
    '''
    data = info[key]
    # union of all adaptors/barcodes (overlapping, but also distinct between samples)
    all_bc = set(bc for sample in data.values() for bc in sample.keys())
    all_bc = sorted(all_bc)

    filename=os.path.join(dir_name, "all_porechop." + key)
    with open(filename, 'w') as file:
        file.write('Sample\t' + '\t'.join(all_bc) + '\n')# print header
        for s in data.keys():
            # list of values for all keys in sample (0 if key does not exist)
            if key=="trimmed":
                sk = list(data[s].keys())
                # Boolean list of trimming status
                sv = [str(k in sk) for k in all_bc]
            else:
                # percentage identities during detection
                sv = [data[s].get(key, 0) for key in all_bc]

            file.write(s + '\t' + '\t'.join(sv) + '\n')


rule porechop_final:
    input: expand_project_path(target_info)
    output: touch("flags/05_porechop.done")
    run:
        # initialize dictionary
        info={}
        info['trimmed']={}
        info['best_start']={}
        info['best_end']={}
        # collect porechop.info for all samples
        for filename in input:
            print('input: ',input)
            get_trimmed(filename, info)
            get_percent(filename, info)
        # write information to file (multiqc compatible)
        # all input files reside in the same directory .../QC/Samples
        # this is the one that will be screened by multiqc
        dir_name=os.path.dirname(input[0])
        if not os.path.exists(dir_name):
            print(f"Error: directory does not exist: {dir_name}")
            exit(1)

        print_info(info, dir_name, key='trimmed')
        print_info(info, dir_name, key="best_start")
        print_info(info, dir_name, key="best_end")

rule qc_porechop:
    input:
        fastq = source
    output:
        info = target_info
    wildcard_constraints:
        # exclude all sample_name that end on "_porechop.info" (already chopped) 
        sample_name = r'(?!.*\.porechop\.info$).*',
    threads: 16
    conda:
        "ont-ppp-porechop"
    params:
        flag = "-v 3",
        nlines = int(config['porechop']['sample_reads']) * 4,
        subset = subset_fastq,
        target = target_fastq
    log:
        logpat
    benchmark:
        bchpat
    shell:'''
        echo "extractting {params.subsample} reads"
        gunzip -c {input.fastq} | head -n {params.nlines} > {params.subset}
        
        ls -lh {params.subset}
        
        echo "porechop_abi {params.flag} -t {threads} -i {params.subset} -o {params.target} > {output.info} 2> {log}"
        porechop_abi {params.flag} -t {threads} -i {params.subset} -o {params.target} > {output.info} 2> {log}
    '''
