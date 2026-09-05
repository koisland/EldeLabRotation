ALN_CFG = config["align"]
ALN_CFG["samples"] = SAMPLES


module Align:
    snakefile:
        "Snakemake-Aligner/workflow/Snakefile"
    config:
        ALN_CFG


use rule * from Align as self_aln_*


rule read_align_all:
    input:
        expand(rules.self_aln_align.input, sm=SAMPLE_NAMES),
