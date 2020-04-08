#!/usr/bin/env python3
import subprocess as sp
def mapping_dna(config):
    """
        mapping DNA using minimap2
    """
    genome_index(config)
    organism = config["data"]["ref"]

    cmd = config["mapping"]["mapping_cmd"]+" "
    cmd += config["mapping"]["mapping_options"]+" "
    cmd += config["data"]["mapping"] + "/" + organism + "_genome.fa "
    cmd += config["info_dict"]["fastq"]+ "/" + config["data"]["Sample_Name"]+".fastq.gz |"
    cmd += config["mapping"]["samtools_cmd"]+" sort "+config["mapping"]["samtools_options"]
    cmd += " -o "+ config["data"]["mapping"] + "/" +config["data"]["Sample_Name"]+".bam ; "
    cmd += config["mapping"]["samtools_cmd"]+ " index "
    cmd += config["data"]["mapping"] + "/" +config["data"]["Sample_Name"]+".bam"
    sp.check_call(cmd, shell=True)

def mapping_rna(config):
    transcriptome_index(config)
    organism = config["data"]["ref"]

    cmd = config["mapping"]["mapping_cmd"]+" "
    cmd += config["mapping"]["mapping_options"]+" "
    cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.fa "
    cmd += config["info_dict"]["fastq"]+ "/" + config["data"]["Sample_Name"]+".fastq.gz |"
    cmd += config["mapping"]["samtools_cmd"]+ " sort "+ config["mapping"]["samtools_options"]
    cmd += " -o "+ config["data"]["mapping"] + "/" +config["data"]["Sample_Name"]+".bam ; "
    cmd += config["mapping"]["samtools_cmd"]+ " index "
    cmd += config["data"]["mapping"] + "/" +config["data"]["Sample_Name"]+".bam"
    print(cmd)
    sp.check_call(cmd, shell=True)

def genome_index(config):
    """
    Generating indices for minimap2
    """
    organism = config["data"]["ref"]
    reference = config["genome"][organism]
    cmd = "ln -s " + reference + " " + config["data"]["mapping"] + "/" + organism + "_genome.fa;"
    cmd += config["mapping"]["mapping_cmd"] +" "+ config["mapping"]["index_options"] + " "
    cmd += config["data"]["mapping"] + "/" + organism + "_genome.mmi "
    cmd += config["data"]["mapping"] + "/" + organism + "_genome.fa "
    sp.check_call(cmd, shell=True)


def transcriptome_index(config):
    """
    Generating indeces for transcriptome
    """
    organism = config["data"]["ref"]
    cmd = config["mapping"]["bedtools_cmd"] +" getfasta "+ config["mapping"]["bedtools_option"]
    cmd += " -fi " + config["genome"][organism]
    cmd += " -bed " + config["transcripts"][organism]
    cmd += " > " + config["data"]["mapping"] + "/" + organism + "_transcripts.fa; "
    cmd += "sed -i \"s/(+)//; s/(-)//\" "
    cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.fa;"
    cmd += config["mapping"]["mapping_cmd"] + " " + config["mapping"]["index_options"] + " "
    cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.mmi "
    cmd += config["data"]["mapping"] + "/" + organism + "_transcripts.fa "
    sp.check_call(cmd, shell=True)
