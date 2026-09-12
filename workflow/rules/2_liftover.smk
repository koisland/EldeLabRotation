REF = config["liftover"]["reference"]["name"]
ANNOTS: dict[str, str] = config["liftover"]["annot"]
MM2_OPTS = config["liftover"]["mm2_opts"]

LIFTOVER_OUTDIR = join(config["output_dir"], "liftover")
LIFTOVER_LOGDIR = join(config["logs_dir"], "liftover")
LIFTOVER_BMKDIR = join(config["benchmarks_dir"], "liftover")

ALN_CFG = {
    "ref": {REF: config["liftover"]["reference"]["path"]},
    "sm": {sm: rules.self_aln_merge_asm_files.output.asm for sm in SAMPLE_NAMES},
    "temp_dir": join(LIFTOVER_OUTDIR, "temp"),
    "output_dir": LIFTOVER_OUTDIR,
    "logs_dir": LIFTOVER_LOGDIR,
    "benchmarks_dir": LIFTOVER_BMKDIR,
    "aln_threads": config["liftover"]["threads_aln"],
    "aln_mem": config["liftover"]["mem_aln"],
    "mm2_opts": MM2_OPTS,
}


wildcard_constraints:
    annot="|".join(ANNOTS.keys()),


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


rule generate_chain_file:
    input:
        expand(rules.asm_ref_bam_to_paf.output, ref=REF, sm="{sm}"),
    output:
        temp(join(LIFTOVER_OUTDIR, REF, "chain", "{sm}.chain")),
    log:
        join(LIFTOVER_LOGDIR, f"chain_{REF}_{{sm}}.log"),
    benchmark:
        join(LIFTOVER_BMKDIR, f"chain_{REF}_{{sm}}.tsv")
    conda:
        "../envs/env.yaml"
    shell:
        """
        paf2chain -i {input} >{output} 2>{log}
        """


rule liftover_genes:
    input:
        gff=lambda wc: ANNOTS[wc.annot],
        target=rules.self_aln_merge_asm_files.output.asm,
        reference=config["liftover"]["reference"]["path"],
    output:
        gff=join(LIFTOVER_OUTDIR, REF, "genes", "{sm}_{annot}.gff"),
    params:
        # https://khchao.com/LiftOn/content/function_manual.html
        min_seq_ident=0.95,
    log:
        join(LIFTOVER_LOGDIR, f"lifton_{REF}_{{sm}}_{{annot}}.log"),
    benchmark:
        join(LIFTOVER_BMKDIR, f"lifton_{REF}_{{sm}}_{{annot}}.tsv")
    conda:
        "../envs/env.yaml"
    shell:
        """
        lifton -g {input.gff} \
        -o {output.gff} \
        -copies \
        -sc {params.min_seq_ident} \
        {input.target} \
        {input.reference}
        """


rule liftover_annotations:
    input:
        # bed or gff file
        annot=lambda wc: ANNOTS[wc.annot],
        chain=rules.generate_chain_file.output,
    output:
        bed=join(LIFTOVER_OUTDIR, REF, "annot", "{sm}_{annot}.bed.gz"),
        unmapped_bed=join(LIFTOVER_OUTDIR, REF, "annot", "{sm}_{annot}_unmapped.bed.gz"),
    log:
        join(LIFTOVER_LOGDIR, f"liftover_{REF}_{{sm}}_{{annot}}.log"),
    benchmark:
        join(LIFTOVER_BMKDIR, f"liftover_{REF}_{{sm}}_{{annot}}.tsv")
    conda:
        "../envs/env.yaml"
    params:
        allow_multiple=(
            "-multiple" if config["liftover"].get("allow_multiple", False) else ""
        ),
        ungzipped_bed=lambda wc, output: output.bed.replace(".gz", ""),
        ungzipped_unmapped_bed=lambda wc, output: output.unmapped_bed.replace(".gz", ""),
    shell:
        """
        liftOver {params.allow_multiple} \
            {input.annot} \
            {input.chain} \
            {params.ungzipped_bed} \
            {params.ungzipped_unmapped_bed} &>{log}
        sort -k1,1 -k2,2n {params.ungzipped_bed} | bgzip >{output.bed}
        sort -k1,1 -k2,2n {params.ungzipped_unmapped_bed} | bgzip >{output.unmapped_bed}
        tabix -p bed {output.bed}
        tabix -p bed {output.unmapped_bed}
        rm -f {params.ungzipped_bed} {params.ungzipped_unmapped_bed}
        """


def liftover_output(sm: str, annot: str):
    infile = ANNOTS[annot]
    if infile.endswith(".gff") or infile.endswith(".gff.gz"):
        return expand(rules.liftover_genes.output, sm=sm, annot=annot)
    elif infile.endswith(".bed"):
        return expand(rules.liftover_annotations.output, sm=sm, annot=annot)
    else:
        raise ValueError(f"Invalid file {infile}. Must be '.bed' or '.gff'.")


rule liftover_all:
    input:
        rules.align_asm_ref_all.input,
        [
            liftover_output(sm, annot)
            for sm in SAMPLE_NAMES
            for annot in ANNOTS.keys()
        ],
