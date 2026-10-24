SEQ_CTX_OUTDIR = join(config["output_dir"], "seq_ctx")
SEQ_CTX_LOGDIR = join(config["logs_dir"], "seq_ctx")
SEQ_CTX_BMKDIR = join(config["benchmarks_dir"], "seq_ctx")


IS_ANNOT_RM = config.get("annot", {}).get("repeatmasker")
OVL_PARAMS = {"prdm9": "-wa -u", "default": "-wa -wb"}
OVL_BEDS = {"prdm9": join(SEQ_CTX_OUTDIR, "{sm}_prdm9_motif.bed")}

if IS_ANNOT_RM:
    OVL_BEDS["repeatmasker"] = join(SEQ_CTX_OUTDIR, "{sm}_rm.bed.gz")


wildcard_constraints:
    annot="|".join(OVL_BEDS.keys()),


# Find PRDM9 motifs
rule find_prdm9_motifs:
    input:
        fa=rules.self_aln_merge_asm_files.output.asm,
    output:
        bed=join(SEQ_CTX_OUTDIR, "{sm}_prdm9_motif.bed"),
    log:
        join(SEQ_CTX_LOGDIR, "find_prdm9_motifs_{sm}.log"),
    benchmark:
        join(SEQ_CTX_BMKDIR, "find_prdm9_motifs_{sm}.tsv")
    conda:
        "../envs/env.yaml"
    params:
        # https://www.nature.com/articles/ng.213
        motif="CCNCCNTNNCCNC",
    shell:
        """
        seqkit locate --ignore-case -d --bed {input.fa} -p {params.motif} >{output.bed} 2>{log}
        """


# RM out file to bed
rule rm_out_to_bed:
    input:
        rm=rules.annot_repeatmasker_output.output if IS_ANNOT_RM else [],
        # https://github.com/logsdon-lab/CenMAP/blob/main/config/repeatmasker_colors.json
        colors=workflow.source_path("../scripts/repeatmasker_colors.json"),
    output:
        rm_bed=join(SEQ_CTX_OUTDIR, "{sm}_rm.bed.gz"),
        rm_bed_idx=join(SEQ_CTX_OUTDIR, "{sm}_rm.bed.gz.tbi"),
    log:
        join(SEQ_CTX_LOGDIR, "overlap_with_rm_{sm}.log"),
    conda:
        "../envs/env.yaml"
    params:
        # https://github.com/logsdon-lab/CenMAP/blob/main/workflow/scripts/create_rm_bed.py
        script=workflow.source_path("../scripts/create_rm_bed.py"),
    shell:
        """
        {{ python {params.script} -i {input.rm} -m {input.colors} \
            | sort -k1,1 -k2,2n \
            | bgzip ;}} >{output.rm_bed} 2>{log}
        tabix -p bed {output.rm_bed} 2>>{log}
        """


# Overlap with liftover annot
rule overlap_with_annot:
    input:
        events=rules.detect_nahr_putative_events.output.bed,
        bed=lambda wc: OVL_BEDS[wc.annot],
    output:
        ovl=join(SEQ_CTX_OUTDIR, "{sm}_{annot}_ovl.bed"),
    log:
        join(SEQ_CTX_LOGDIR, "overlap_with_{sm}_{annot}.log"),
    conda:
        "../envs/env.yaml"
    params:
        ovl_params=lambda wc: OVL_PARAMS.get(wc.annot, OVL_PARAMS["default"]),
    shell:
        """
        bedtools intersect -a {input.events} -b {input.bed} {params.ovl_params} >{output.ovl} 2>>{log}
        """


# Check for microhomology (>=5) at breakpoints per region
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5237343/
rule msa_breakpoint_boundaries:
    input:
        fa=rules.self_aln_merge_asm_files.output.asm,
        events=rules.detect_nahr_putative_events.output.bed,
    output:
        events_scored=join(SEQ_CTX_OUTDIR, "{sm}_w_homology_score.bed"),
    log:
        join(SEQ_CTX_LOGDIR, "msa_breakpoint_boundaries_{sm}.log"),
    benchmark:
        join(SEQ_CTX_BMKDIR, "msa_breakpoint_boundaries_{sm}.tsv")
    conda:
        "../envs/env.yaml"
    params:
        bp_slop=5,
        script=workflow.source_path("../scripts/msa_breakpoint.py"),
    shell:
        """
        python {params.script} -i {input.events} -f {input.fa} -s {params.bp_slop} >{output.events_scored} 2>{log}
        echo "Matches >= 5" >>{log}
        awk '$5>=5' {output.events_scored} | wc -l >>{log}
        echo "Total" >>{log}
        wc -l {input.events} >>{log}
        """


rule sequence_context_all:
    input:
        expand(rules.find_prdm9_motifs.output, sm=SAMPLE_NAMES),
        expand(rules.msa_breakpoint_boundaries.output, sm=SAMPLE_NAMES),
        expand(rules.rm_out_to_bed.output, sm=SAMPLE_NAMES) if IS_ANNOT_RM else [],
        expand(rules.overlap_with_annot.output, sm=SAMPLE_NAMES, annot=OVL_BEDS.keys()),
