ALIGN_OUTDIR = join(config["output_dir"], "align")
ALIGN_LOGDIR = join(config["logs_dir"], "align")
ALIGN_BMKDIR = join(config["benchmarks_dir"], "align")


align_outputs = []
if "align" in config:
    ALN_CFG = config["align"]
    ALN_CFG["samples"] = SAMPLES
    ALN_CFG["output_dir"] = ALIGN_OUTDIR
    ALN_CFG["logs_dir"] = ALIGN_LOGDIR
    ALN_CFG["benchmarks_dir"] = ALIGN_BMKDIR

    module Align:
        snakefile:
            "Snakemake-Aligner/workflow/Snakefile"
        config:
            ALN_CFG

    use rule * from Align as self_aln_*

    align_outputs.extend(expand(rules.self_aln_align.input, sm=SAMPLE_NAMES))

else:

    rule self_aln_merge_asm_files:
        input:
            asm=lambda wc: SAMPLE_INFO[wc.sm]["asm_fa"],
        output:
            asm=join(ALIGN_OUTDIR, "{sm}.fa"),
            idx=join(ALIGN_OUTDIR, "{sm}.fa.fai"),
        conda:
            "../envs/env.yaml"
        shell:
            """
            zcat -f {input.asm} > {output.asm}
            samtools faidx {output.asm}
            """

    rule self_aln_merge_read_asm_alignments:
        input:
            alignment=lambda wc: SAMPLE_INFO[wc.sm]["bam"],
            alignment_idx=lambda wc: SAMPLE_INFO[wc.sm]["bam"] + ".bai",
        output:
            alignment=join(ALIGN_OUTDIR, "{sm}.bam"),
            alignment_idx=join(ALIGN_OUTDIR, "{sm}.bam.bai"),
        shell:
            """
            ln -s {input.alignment} {output.alignment}
            ln -s {input.alignment_idx} {output.alignment_idx}
            """

    align_outputs.extend(
        expand(rules.self_aln_merge_asm_files.output, sm=SAMPLE_NAMES)
    )
    align_outputs.extend(
        expand(rules.self_aln_merge_read_asm_alignments.output, sm=SAMPLE_NAMES)
    )


rule read_align_all:
    input:
        align_outputs,
