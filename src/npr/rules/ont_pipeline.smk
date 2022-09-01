#snakefile
from npr.snakehelper import retRule

include: retRule("1_basecalling.smk", config)
include: retRule("2_rename.smk", config)
include: retRule("3_qc.smk", config)
include: retRule("4_data_transfer.smk", config)
include: retRule("5_mapping.smk", config)

rule all:
    input:
        'flags/1_basecalling.done',
        'flags/2_renamed.done',
        'flags/3_qc.done',
        'flags/4_transfer.done',
        'flags/5_mapping.done'
