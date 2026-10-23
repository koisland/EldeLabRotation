
SEQ_CTX_OUTDIR = join(config["output_dir"], "seq_ctx")
SEQ_CTX_LOGDIR = join(config["logs_dir"], "seq_ctx")
SEQ_CTX_BMKDIR = join(config["benchmarks_dir"], "seq_ctx")


# Find and overlap PRDM9 motifs
rule find_overlap_prdm9_motifs:
    input:
        fa=rules.self_aln_merge_asm_files.output.asm,
        events=rules.detect_nahr_putative_events.output.bed,
    output:
        bed=join(SEQ_CTX_OUTDIR, "{sm}_prdm9_motif.bed"),
        ovl=join(SEQ_CTX_OUTDIR, "{sm}_prdm9_motif_ovl.bed"),
    log:
        join(SEQ_CTX_LOGDIR, "find_overlap_prdm9_motifs_{sm}.log"),
    benchmark:
        join(SEQ_CTX_BMKDIR, "find_overlap_prdm9_motifs_{sm}.tsv")
    params:
        # https://www.nature.com/articles/ng.213
        motif="CCNCCNTNNCCNC",
    conda:
        "../envs/env.yaml",
    shell:
        """
        seqkit locate --ignore-case -d --bed {input.fa} -p {params.motif} >{output.bed} 2>{log}
        bedtools intersect -a {input.events} -b {output.bed} -wa -u >{output.ovl} 2>>{log}
        """


# Overlap with RM
rule overlap_with_rm:
    input:
        rm="",
        events=rules.detect_nahr_putative_events.output.bed,
    output:
        bed="",
    conda:
        "../envs/env.yaml",
    shell:
        """
        """

# Overlap with liftover annot

rule msa_breakpoint_boundaries:
    input:
        fa=rules.self_aln_merge_asm_files.output.asm,
        events=rules.detect_nahr_putative_events.output.bed,
    output:
        msa=""
    conda:
        "../envs/env.yaml",
    shell:
        """
        mafft
        """


rule sequence_context_all:
    input:
        expand(rules.find_overlap_prdm9_motifs.output, sm=SAMPLE_NAMES)
