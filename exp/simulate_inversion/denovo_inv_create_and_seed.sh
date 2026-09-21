#!/bin/bash

set -euo pipefail

sample="${1}"

# Simulate inversion
fa_out_inv_hap="exp/simulate_inversion/${sample}_hap1_after_inv.fa"
python exp/simulate_inversion/denovo_inv_create.py \
    -f "exp/simulate_inversion/${sample}_TCAF_hap1.fa" \
    -s 42 \
    -o "exp/simulate_inversion/${sample}_hap1"

# Generate simulated hifi reads from inversion
# Make ideal as possible
fq_out_inv_sim="exp/simulate_inversion/${sample}_hap1_sim.fq.gz"
badread simulate \
    --reference "${fa_out_inv_hap}" \
    --quantity 5x \
    --error_model random \
    --qscore_model ideal \
    --glitches 0,0,0 \
    --junk_reads 0 \
    --random_reads 0 \
    --chimeras 0 \
    --identity 30,3 \
    --length 20000,2000 \
    --start_adapter_seq "" \
    --end_adapter_seq "" \
    --seed 42 \
    | bgzip > "${fq_out_inv_sim}"

# Remove badread comments
# Add SM: tag to pick out
seqkit replace -p "\s.+" "${fq_out_inv_sim}" \
    | seqkit replace -p $ -r " SM:Z:inv_${sample}" \
    | bgzip > "${fq_out_inv_sim}.tmp"
mv "${fq_out_inv_sim}.tmp" "${fq_out_inv_sim}"

# Align simulated inv hap to region
bash exp/simulate_inversion/existing_haps_align_rgn.sh \
    "results/align/${sample}.bam" \
    "results/align/${sample}.fa" \
    "exp/simulate_inversion/${sample}_TCAF_regions.bed" \
    "${fq_out_inv_sim}"
