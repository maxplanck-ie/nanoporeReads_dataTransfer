import os
import sys
import subprocess as sp
import warnings

def call_cmd(command):
    sp.check_output(command, shell=True)

rule bc_split:
    input:
        "demux.done"
    output:
        touch("bc.split")
    run:
        if config["bc_kit"] != "no_bc":
            cmd = config["pycoQc"]["barcodeSplit"]
            cmd += "fastq/sequencing_summary.txt "
            cmd += "-o fastq/"
            sp.check_call(cmd, shell=True)

rule file_rename:
    input:
        bc_split = "bc.split",
    output:
        rename = temp(touch("{sample_id}_renamed.done"))
    log:
        out ="LOG/rename_{sample_id}.log.out",
        err = "LOG/rename_{sample_id}.log.err"
    run:
        this_sample = config["data"][wildcards.sample_id]
        if this_sample["barcode_kits"] != "no_bc":
            bs = this_sample["index_id"]
            bs = [s for s in " ".join(list(bs)).split() if s.isdigit()]
            bs = "".join(bs)
            print("bs:", bs)
            bs_fastq = os.path.join("fastq","barcode"+bs)
            this_sample["index_id"] = "barcode"+bs
        else:
            assert(len(config["data"])==1)
            bs_fastq = os.path.join("fastq")

        if not os.path.exists("Project_"+this_sample["Sample_Project"]):
           os.mkdir("Project_"+this_sample["Sample_Project"])
        if not os.path.exists("Project_"+this_sample["Sample_Project"]+"/Sample_"+this_sample["Sample_ID"]):
            os.mkdir("Project_"+this_sample["Sample_Project"]+"/Sample_"+this_sample["Sample_ID"])
            sample_path = os.path.join("Project_"+this_sample["Sample_Project"]+"/Sample_"+this_sample["Sample_ID"])
            sample_name = this_sample["Sample_Name"]
            try:
                os.path.exists(os.path.join("fastq", "pass"))
                bs_fastq_pass = os.path.join("fastq", "pass")
                bs_fastq_fail = os.path.join("fastq", "fail")
                if this_sample["barcode_kits"] != "no_bc":
                    bs_fastq_pass = os.path.join("fastq", "pass","barcode"+bs)
                    bs_fastq_fail = os.path.join("fastq", "fail","barcode"+bs)
                cmd = "cat {}/*.fastq.gz {}/*.fastq.gz > {}/{}.fastq.gz;".format(bs_fastq_pass, bs_fastq_fail,sample_path,sample_name)
            except:
                cmd = "cat {}/*.fastq.gz > {}/{}.fastq.gz;".format(bs_fastq,sample_path,sample_name)

            cmd += " >> "+log.out+" 2> "+log.err
            with open(log.out, "w") as called_cmd:
                called_cmd.write(cmd)
            called_cmd.close()
            sp.check_call(cmd, shell=True)
        else:
            warnings.warn("Directory called: Project_"+this_sample["Sample_Project"]+\
                        "/Sample_"+this_sample["Sample_ID"]+" already exists. It stays untouched!")
        #TODO Do we want to keep or remove all the fastq files with wrong or unknown barcodes?


rule rename_done:
    input:
        expand("{sample_id}_renamed.done", sample_id = config["data"].keys())
    output:
        touch("renamed.done")
