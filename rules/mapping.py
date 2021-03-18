#!/usr/bin/env python3
import subprocess as sp
import os

def genome_index(config, path):
    """
    Generating genome indices for minimap2
    """
    reference = config["organism"]
    cmd = "ln -s " + reference + " " + path + "/" + reference + "_genome.fa;"
    cmd += config["mapping"]["mapping_cmd"] +" "+ config["mapping"]["index_options"] + " "
    cmd += path + "/" + reference + "_genome.mmi "
    cmd += path + "/" + reference + "_genome.fa "
    print(cmd)
    sp.check_call(cmd, shell=True)




def mapping_rna_contamination(config, data):
   """
   Mapping RNA/cDNA using minimap2 to several genomes and rRNA to look for hte contamiantion
   """
   # TODO different organism on the same flowcell or different projects!
   h_transcripts = config["transcripts"]["hg38"]
   m_transcripts = config["transcripts"]["mm10"]
   os.mkdir(config["info_dict"]["flowcell_path"]+"/contamination_report")
   report_dir = config["info_dict"]["flowcell_path"]+"/contamination_report/"
   for ref in ["hg38","mm10", "human_rRNA", "mouse_rRNA"]:
        genome_index(config,ref, report_dir)
   for k, v in data.items():
        bams = ""
        names = ""
        analysis_dir = report_dir+v["Sample_Project"]
        if not os.path.exists(analysis_dir):
            os.mkdir(analysis_dir)

        if not os.path.exists(analysis_dir+"/"+v["Sample_ID"]):
            os.mkdir(analysis_dir+"/"+v["Sample_ID"])
        for index, ref in enumerate(["hg38","mm10"]):
            transcripts = h_transcripts
            organism = "hGenome"
            if index == 1:
                transcripts= m_transcripts
                organism = "mGenome"
            cmd = config["mapping"]["mapping_cmd"]+ " "
            cmd += config["mapping"]["mapping_rna_options"]
            cmd += " --junc-bed "+transcripts + " "
            cmd += report_dir + ref + "_genome.fa "
            cmd += config["info_dict"]["flowcell_path"]+"/Project_"+v["Sample_Project"]+"/Sample_"+v["Sample_ID"]+"/"+v["Sample_Name"]+".fastq.gz | "
            cmd += config["mapping"]["samtools_cmd"]+ " sort "+ config["mapping"]["samtools_options"]
            cmd += " -o "+ analysis_dir + "/"+v["Sample_ID"]+"/" +v["Sample_Name"]+"."+organism+".bam ; "
            cmd += config["mapping"]["samtools_cmd"]+ " index "
            cmd += analysis_dir + "/"+v["Sample_ID"]+"/" +v["Sample_Name"]+"."+organism+".bam"
            print(cmd)
            bams += analysis_dir + "/"+v["Sample_ID"]+"/" +v["Sample_Name"]+"."+organism+".bam "
            names += organism+" "
            sp.check_call(cmd, shell=True)
        for index, ref in enumerate(["human_rRNA", "mouse_rRNA"]):
            organism = "hrRNA"
            if index == 1:
                organism = "mrRNA"
            cmd = config["mapping"]["mapping_cmd"]+" "
            cmd += config["mapping"]["mapping_dna_options"]+" "
            cmd += report_dir + ref + "_genome.fa "
            cmd += config["info_dict"]["flowcell_path"]+"/Project_"+v["Sample_Project"]+"/Sample_"+v["Sample_ID"]+"/"+v["Sample_Name"]+".fastq.gz | "
            cmd += config["mapping"]["samtools_cmd"]+ " sort "+ config["mapping"]["samtools_options"]
            cmd += " -o "+ analysis_dir + "/"+v["Sample_ID"]+"/" +v["Sample_Name"]+"."+organism+".bam ; "
            cmd += config["mapping"]["samtools_cmd"]+ " index "
            cmd += analysis_dir + "/"+v["Sample_ID"]+"/" +v["Sample_Name"]+"."+organism+".bam"
            print(cmd)
            bams += analysis_dir + "/"+v["Sample_ID"]+"/" +v["Sample_Name"]+"."+organism+".bam "
            names += organism+" "
            sp.check_call(cmd, shell=True)

        cmd = config["nanocomp"]["qc_cmd"]
        cmd += " --bam {}".format(bams)
        cmd += " -o {} ".format(analysis_dir + "/"+v["Sample_ID"]+"/")
        cmd += config["nanocomp"]["qc_options"]+" "+names
        print(cmd)
        sp.check_call(cmd, shell=True)



# def transcriptome_index(config):
#     """
#     Generating indeces for transcriptome
#     """
#     organism = config["data"]["ref"]
#     cmd = config["mapping"]["bedtools_cmd"] +" getfasta "+ config["mapping"]["bedtools_option"]
#     cmd += " -fi " + config["genome"][organism]
#     cmd += " -bed " + config["transcripts"][organism]
#     cmd += " > " + config["data"]["mapping"] + "/" + organism + "_transcripts.fa; "
#     cmd += "sed -i \"s/(+)//; s/(-)//\" "
#     cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.fa;"
#     cmd += config["mapping"]["mapping_cmd"] + " " + config["mapping"]["index_options"] + " "
#     cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.mmi "
#     cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.fa "
#     sp.check_call(cmd, shell=True)
#
# def mapping_rna(config):
#     """
#     mapping RNA using minimap2
#     """
#     transcriptome_index(config)
#     organism = config["data"]["ref"]
#     cmd = config["mapping"]["mapping_cmd"]+" "
#     cmd += config["mapping"]["mapping_options"]+" "
#     cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.fa "
#     cmd += config["info_dict"]["fastq"]+ "/" + config["data"]["Sample_Name"]+".fastq.gz |"
#     cmd += config["mapping"]["samtools_cmd"]+ " sort "+ config["mapping"]["samtools_options"]
#     cmd += " -o "+ config["data"]["mapping"] + "/" +config["data"]["Sample_Name"]+".bam ; "
#     cmd += config["mapping"]["samtools_cmd"]+ " index "
#     cmd += config["data"]["mapping"] + "/" +config["data"]["Sample_Name"]+".bam"
#     print(cmd)
#     sp.check_call(cmd, shell=True)


rule mapping_data:
    input:
        transferred = "{sample_id}.transferred"
    output:
        mapped = "{sample_id}.bam"
    run:
        group=this_sample["Sample_Project"].split("_")[-1].lower()
        final_path = os.path.join(config["paths"]["groupDir"],group,"sequencing_data/OxfordNanopore/"+config["input"]["name"])
        analysis = os.path.join(final_path,"Analysis_"+this_sample["Sample_Project"])
        analysis = analysis+"/mapping_on_"+config["organism"]
        genome_index(config, analysis)
        mapping_path = os.path.join(analysis, "/Sample_"+wildcards.sample_id)
        cmd = config["mapping"]["mapping_cmd"]+ " "
        if ("RNA" in config["info_dict"]["kit"]) or (config["protocol"] == 'rna') or (config["protocol"] == 'cdna'):
            cmd += config["mapping"]["mapping_rna_options"]
            cmd += " --junc-bed "+config["transcripts"][config["organism"]] + " "
        else:
            cmd += config["mapping"]["mapping_dna_options"]+" "
        cmd += analysis_dir + "/" + config["organism"] + "_genome.fa "
        cmd += config["info_dict"]["flowcell_path"]+"/Project_"+this_sample["Sample_Project"]+"/Sample_"+this_sample["Sample_ID"]+"/"+this_sample["Sample_Name"]+".fastq.gz | "
        cmd += config["mapping"]["samtools_cmd"]+ " sort "+ config["mapping"]["samtools_options"]
        cmd += " -o "+ mapping_path+"/" +this_sample["Sample_Name"]+".bam ; "
        cmd += config["mapping"]["samtools_cmd"]+ " index "
        cmd += mapping_path+"/" +this_sample["Sample_Name"]+".bam"
        print(cmd)
        sp.check_call(cmd, shell=True)
