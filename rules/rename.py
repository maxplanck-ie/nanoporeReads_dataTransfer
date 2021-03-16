import os
import sys
import subprocess as sp
import warnings


rule bc_split:
    input:
        fastq = "fastq"
    output:
        file = temp("bc.split")
    run:
        if config["bc_kit"] != "no_bc":
            cmd = config["pycoQc"]["barcodeSplit"]
            cmd += "fastq/sequencing_summary.txt "
            cmd += "-o fastq/; touch {output.file}"
            print(cmd)
            sp.check_call(cmd, shell=True)
        else:
            sp.check_call("touch {output.file}", shell=True)


rule file_rename:
    input:
        bc_split = "bc.split",
        fastq = "fastq"
    output:
        rename = temp("{sample_id}_renamed.done")
    run:
        this_sample = config["data"][wildcards.sample_id]
        if this_sample["barcode_kits"] is not "no_bc":
            bs = this_sample["index_id"]
            bs = [s for s in " ".join(list(bs)).split() if s.isdigit()]
            bs = "".join(bs)
            print("bs:", bs)
            bs_fastq = os.path.join(input.fastq,"barcode"+bs)
            this_sample["index_id"] = "barcode"+bs
        else:
            assert(len(config["data"])==1)

        if not os.path.exists("Project_"+this_sample["Sample_Project"]):
           os.mkdir("Project_"+this_sample["Sample_Project"])
        if not os.path.exists("Project_"+this_sample["Sample_Project"]+"/Sample_"+this_sample["Sample_ID"]):
            os.mkdir("Project_"+this_sample["Sample_Project"]+"/Sample_"+this_sample["Sample_ID"])
            sample_path = os.path.join("Project_"+this_sample["Sample_Project"]+"/Sample_"+this_sample["Sample_ID"])
            sample_name = this_sample["Sample_Name"]
            cmd = "cat {}/*.fastq.gz > {}/{}.fastq.gz ;".format(bs_fastq,sample_path,sample_name)
            sp.check_call(cmd, shell=True)
        else:
            warnings.warn("Directory called: Project_"+this_sample["Sample_Project"]+\
                        "/Sample_"+this_sample["Sample_ID"]+" already exists. It stays untouched!")
        #TODO Do we want to keep or remove all the fastq files with wrong or unknown barcodes?
        cmd = "touch "+output.rename
        sp.check_call(cmd, shell=True)
