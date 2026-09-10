ANNOTS = config["liftover"]["annot"].keys()
MM2_OPTS = config["liftover"]["mm2_opts"]

LIFTOVER_OUTDIR = config["liftover"]["output_dir"]
LIFTOVER_LOGDIR = config["liftover"]["logs_dir"]
LIFTOVER_BMKDIR = config["liftover"]["benchmarks_dir"]

ALN_CFG = {
    "ref": {
        config["liftover"]["reference"]["name"]: config["liftover"]["reference"]["path"]
    },
    "sm": {sm: rules.self_aln_merge_asm_files.output.asm for sm in SAMPLE_NAMES},
    "temp_dir": join(LIFTOVER_OUTDIR, "temp"),
    "output_dir": LIFTOVER_OUTDIR,
    "logs_dir": LIFTOVER_LOGDIR,
    "benchmarks_dir": LIFTOVER_BMKDIR,
    "aln_threads": config["liftover"]["threads_aln"],
    "aln_mem": config["liftover"]["mem_aln"],
    "mm2_opts": MM2_OPTS,
}


module AlignToRef:
    snakefile:
        "asm-to-reference-alignment/workflow/Snakefile"
    config:
        ALN_CFG


use rule * from AlignToRef as asm_ref_*


use rule all from AlignToRef as align_asm_ref_all with:
    default_target: True
    input:
        rules.asm_ref_all.input,


# TODO: Generate chain file

# TODO: Liftover


rule liftover_all:
    input:
        rules.align_asm_ref_all.input,
