ALN_CFG = config["align"]


module Align:
    snakefile:
        github(
            "logsdon-lab/Snakemake-Aligner",
            path="workflow/Snakefile",
            commit="18b558fe551b12186ee74b535fb37849d3d36572",
        )
    config:
        ALN_CFG


use rule * from Align as self_aln_*
