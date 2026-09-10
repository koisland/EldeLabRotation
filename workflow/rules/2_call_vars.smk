CALL_VAR_OUTDIR = config["call_variants"]["output_dir"]
CALL_VAR_LOGDIR = config["call_variants"]["logs_dir"]
CALL_VAR_BMKDIR = config["call_variants"]["benchmarks_dir"]


# Sniffles mosaic
rule call_variants:
    input:
        bam=rules.self_aln_merge_read_asm_alignments.output.alignment,
        ref=rules.self_aln_merge_asm_files.output.asm,
    output:
        vcf=join(CALL_VAR_OUTDIR, "{sm}.vcf"),
    conda:
        "../envs/env.yaml"
    threads: config["call_variants"]["threads_sniffles"]
    resources:
        mem=config["call_variants"]["mem_sniffles"],
    shell:
        """
        sniffles \
            --input {input.bam} \
            --reference {input.ref} \
            --vcf {output.vcf} \
            --mosaic \
            --threads {threads} \
            --output-rnames
        """


# HiPhase
rule phase_variants_bam:
    input:
        bam=rules.self_aln_merge_read_asm_alignments.output.alignment,
        vcf=rules.call_variants.output.vcf,
        ref=rules.self_aln_merge_asm_files.output.asm,
    output:
        hvcf=join(CALL_VAR_OUTDIR, "{sm}.phased.vcf"),
        hbam=join(CALL_VAR_OUTDIR, "{sm}.phased.bam"),
    conda:
        "../envs/env.yaml"
    threads: config["call_variants"]["threads_hiphase"]
    resources:
        mem=config["call_variants"]["mem_hiphase"],
    shell:
        """
        hiphase \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --output-vcf {output.hvcf} \
            --output-bam {output.hbam} \
            --reference {input.ref} \
            --threads {threads}
        """


rule call_vars_all:
    input:
        expand(rules.call_variants.output, sm=SAMPLE_NAMES),
        expand(rules.phase_variants_bam.output, sm=SAMPLE_NAMES),
