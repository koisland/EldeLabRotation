ANNOT_OUTDIR = config["annot"]["output_dir"]
SPLIT_MULTIFA_DIR = join(config["annot"]["output_dir"], "fa")
ANNOT_LOGDIR = config["annot"]["logs_dir"]
ANNOT_BMKDIR = config["annot"]["benchmarks_dir"]

ANNOT_CFG = {
    "samples": {
        sm: {"fa": rules.self_aln_merge_asm_files.output.asm} for sm in SAMPLE_NAMES
    },
    "output_dir": ANNOT_OUTDIR,
    "log_dir": ANNOT_LOGDIR,
    "benchmark_dir": ANNOT_BMKDIR,
    "repeatmasker": {
        "species": config["annot"]["repeatmasker_species"],
        "engine": config["annot"]["repeatmasker_engine"],
        "threads": config["annot"]["threads_annot"],
        "mem": config["annot"]["mem_annot"],
    },
}


module AnnotateRepeats:
    snakefile:
        "Snakemake-Repeat-Annotation/workflow/Snakefile"
    config:
        ANNOT_CFG


use rule * from AnnotateRepeats as annot_*


# rule run_biser:
#     input:
#         ""
#     output:
#         ""
#     conda:
#         "../envs/env.yaml"
#     shell:
#         """
#         """


rule annote_ampl_all:
    input:
        rules.annot_all.input,
