ANNOT_OUTDIR = join(config["output_dir"], "annot")
SPLIT_MULTIFA_DIR = join(ANNOT_OUTDIR, "fa")
ANNOT_LOGDIR = join(config["logs_dir"], "annot")
ANNOT_BMKDIR = join(config["benchmarks_dir"], "annot")


ANNOT_CFG = {
    "samples": {
        sm: {"fa": rules.self_aln_merge_asm_files.output.asm} for sm in SAMPLE_NAMES
    },
    "output_dir": ANNOT_OUTDIR,
    "log_dir": ANNOT_LOGDIR,
    "benchmark_dir": ANNOT_BMKDIR,
    "repeatmasker": {
        "species": config["annot"]["species_rm"],
        "engine": config["annot"]["engine_rm"],
        "threads": config["annot"]["threads_rm"],
        "mem": config["annot"]["mem_rm"],
    },
    "biser": {
        "threads": config["annot"]["threads_biser"],
        "mem": config["annot"]["mem_biser"],
    },
}


# TODO: TRF + windowmasker to better mask reference. Look at Eichler Lab GH
# https://github.com/EichlerLab/sedef_smk/tree/main/rules
# https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=4161691683_uqKBa5fu135H1RWEk6zfhhtEF4WH&db=hub_4837794_T2T-CHM13v2.0&c=chr7&g=hub_4837794_segDups2024
module AnnotateRepeats:
    snakefile:
        "Snakemake-Repeat-Annotation/workflow/Snakefile"
    config:
        ANNOT_CFG


use rule * from AnnotateRepeats as annot_*


rule annote_ampl_all:
    input:
        rules.annot_all.input,
