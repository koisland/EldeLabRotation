#!/bin/bash

set -euo pipefail

fa=$1
outdir=$2

split_dir="${outdir}/fasta"
sim_dir="${outdir}/sim"
rm -rf "${split_dir}" "${sim_dir}"
mkdir -p "${sim_dir}"

seqkit split -i -O "${split_dir}" "${fa}" --by-id-prefix "split_"
for file in "${split_dir}"/split_*.gz; do
  bname=$(basename "${file_nogz}")
  badread simulate \
    --reference "${file}" \
    --quantity 1x \
    --error_model random \
    --qscore_model ideal \
    --glitches 0,0,0 \
    --junk_reads 0 \
    --random_reads 0 \
    --chimeras 0 \
    --length 20000,2000 \
    --start_adapter_seq "" \
    --end_adapter_seq "" \
    --identity 30,3 | bgzip > "${sim_dir}/${bname}"
done
