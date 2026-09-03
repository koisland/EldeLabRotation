ALN_CFG = config["align"]


module Align:
    snakefile:
        "Snakemake-Aligner/workflow/Snakefile"
    config:
        ALN_CFG


use rule * from Align as self_aln_*
