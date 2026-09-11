REF = config["liftover"]["reference"]["name"]
ANNOTS: dict[str, str] = config["liftover"]["annot"]
MM2_OPTS = config["liftover"]["mm2_opts"]

LIFTOVER_OUTDIR = config["liftover"]["output_dir"]
LIFTOVER_LOGDIR = config["liftover"]["logs_dir"]
LIFTOVER_BMKDIR = config["liftover"]["benchmarks_dir"]

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
        paf2chain -i {input} >{output} &>{log}
        """


rule liftover_annotations:
    input:
        # bigbed file
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
        """


rule liftover_all:
    input:
        rules.align_asm_ref_all.input,
        expand(rules.liftover_annotations.output, sm=SAMPLE_NAMES, annot=ANNOTS.keys()),
