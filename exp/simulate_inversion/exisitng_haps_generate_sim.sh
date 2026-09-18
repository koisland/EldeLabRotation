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
    # https://github.com/yukiteruono/pbsim3/blob/master/data/ERRHMM-RSII.model
    gunzip ${file} || true
    file_nogz=$(echo "${file}" | sed 's/.gz//g')
    bname_file_nogz=$(basename "${file_nogz}" .fa | sed -e 's/split_//g' -e 's/__/:/g')
    pbsim --strategy wgs \
      --method qshmm \
      --qshmm exp/simulate_inversion/QSHMM-RSII.model \
      --depth 5 \
      --genome "${file_nogz}"  \
      --prefix "${sim_dir}/${bname_file_nogz}"
done
