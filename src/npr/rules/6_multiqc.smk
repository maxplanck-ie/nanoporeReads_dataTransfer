#!/usr/bin/env python3
import os
import warnings


rule final_multiqc:
    input:
        "flags/5_mapping.done",
    output:
        touch("flags/6_multiqc.done"),
    params:
        finalpath=config["data"]["finalpath"],
        multiqc_dir=os.path.join(finalpath, "multiqc"),
    log:
        err="log/multiqc/project-final_multiqc.err",
    shell:
        """
        mkdir -p {params.multiqc_dir}
        multiqc -o {params.multiqc_dir} {params.finalpath} 2> {log.err}
    """
