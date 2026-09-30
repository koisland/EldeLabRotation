CALL_NAHR_OUTDIR = join(config["output_dir"], "call_nahr")
CALL_NAHR_LOGDIR = join(config["logs_dir"], "call_nahr")
CALL_NAHR_LOGDIR = join(config["benchmarks_dir"], "call_nahr")


# workflow/scripts/RECallS/target/release/RECallS
rule detect_nahr_putative_events:
    input:
        bam=rules.self_aln_merge_read_asm_alignments.output.alignment,
        ref=rules.self_aln_merge_asm_files.output.asm,
        bn=config["call_nahr"]["recalls_bin"],
        ignore_bed=(rules.annot_longdust.output.bed if config.get("annot") else []),
    output:
        output_inv_bed=join(CALL_NAHR_OUTDIR, "{sm}", "calls_inv.bed"),
        output_del_bed=join(CALL_NAHR_OUTDIR, "{sm}", "calls_del.bed"),
    log:
        join(CALL_NAHR_LOGDIR, "detect_nahr_put_events_{sm}.log"),
    benchmark:
        join(CALL_NAHR_LOGDIR, "detect_nahr_put_events_{sm}.tsv")
    conda:
        "../envs/env.yaml"
    threads: config["call_nahr"]["threads"]
    resources:
        mem=config["call_nahr"]["mem"],
    params:
        output_dir=lambda wc, output: dirname(str(output.output_inv_bed)),
        ignore_bed=lambda wc, input: (
            f"-n {input.ignore_bed}" if input.ignore_bed else ""
        ),
    shell:
        """
        ./{input.bn} \
            -i {input.bam} \
            -f {input.ref} \
            -o {params.output_dir} \
            -t {threads} \
            {params.ignore_bed} 2>{log}
        """


rule detect_nahr_events_all:
    input:
        expand(rules.detect_nahr_putative_events.output, sm=SAMPLE_NAMES),
