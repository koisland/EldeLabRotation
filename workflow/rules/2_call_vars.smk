CALL_VAR_OUTDIR = config["call_variants"]["output_dir"]
CALL_VAR_LOGDIR = config["call_variants"]["logs_dir"]
CALL_VAR_BMKDIR = config["call_variants"]["benchmarks_dir"]


rule run_deepvariant:
    input:
        bam=rules.self_aln_merge_read_asm_alignments.output.alignment,
        ref=rules.self_aln_merge_asm_files.output.asm,
    output:
        vcf=join(CALL_VAR_OUTDIR, "{sm}_deepvariant.gvcf.gz"),
        # https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
        gvcf=join(CALL_VAR_OUTDIR, "{sm}_deepvariant.vcf.gz"),
    log:
        join(CALL_VAR_LOGDIR, "run_deepvariant_{sm}.log"),
    benchmark:
        join(CALL_VAR_BMKDIR, "run_deepvariant_{sm}.tsv")
    singularity:
        "docker://google/deepvariant:1.9.0"
    threads: config["call_variants"]["threads_deepvariant"]
    resources:
        mem=config["call_variants"]["mem_deepvariant"],
    params:
        model="PACBIO",
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
            --model_type {params.model} \
            --ref {input.ref} \
            --reads {input.bam} \
            --output_vcf {output.vcf} \
            --output_gvcf {output.gvcf} \
            --num_shards {threads} &>{log}
        """


# Sniffles mosaic
rule run_sniffles:
    input:
        bam=rules.self_aln_merge_read_asm_alignments.output.alignment,
        ref=rules.self_aln_merge_asm_files.output.asm,
    output:
        vcf=join(CALL_VAR_OUTDIR, "{sm}_sniffles.vcf.gz"),
    log:
        join(CALL_VAR_LOGDIR, "run_sniffles_{sm}.log"),
    benchmark:
        join(CALL_VAR_BMKDIR, "run_sniffles_{sm}.tsv")
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
            --output-rnames &>{log}
        """


# HiPhase
rule phase_variants_bam:
    input:
        bam=rules.self_aln_merge_read_asm_alignments.output.alignment,
        deepvariant_vcf=rules.run_deepvariant.output.vcf,
        sniffles_vcf=rules.run_sniffles.output.vcf,
        ref=rules.self_aln_merge_asm_files.output.asm,
    output:
        deepvariant_hvcf=join(CALL_VAR_OUTDIR, "{sm}_deepvariant.phased.vcf.gz"),
        sniffles_hvcf=join(CALL_VAR_OUTDIR, "{sm}_sniffles.phased.vcf.gz"),
        hbam=join(CALL_VAR_OUTDIR, "{sm}.phased.bam"),
    log:
        join(CALL_VAR_LOGDIR, "hiphase_{sm}.log"),
    benchmark:
        join(CALL_VAR_BMKDIR, "hiphase_{sm}.tsv")
    conda:
        "../envs/env.yaml"
    threads: config["call_variants"]["threads_hiphase"]
    resources:
        mem=config["call_variants"]["mem_hiphase"],
    shell:
        """
        hiphase \
            --bam {input.bam} \
            --vcf {input.deepvariant_vcf} \
            --output-vcf {output.deepvariant_hvcf} \
            --vcf {input.sniffles_vcf} \
            --output-vcf {output.sniffles_hvcf} \
            --output-bam {output.hbam} \
            --reference {input.ref} \
            --threads {threads} \
            --ignore-read-groups &>{log}
        """


rule call_vars_all:
    input:
        expand(rules.run_sniffles.output, sm=SAMPLE_NAMES),
        expand(rules.run_deepvariant.output, sm=SAMPLE_NAMES),
        expand(rules.phase_variants_bam.output, sm=SAMPLE_NAMES),
