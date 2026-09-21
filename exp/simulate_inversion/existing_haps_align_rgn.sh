#!/bin/bash

set -euxo pipefail

wd=$(dirname "${0}")
bam="${1}"
asm="${2}"
regions_file="${3}"
add_fa="${4}"

reads="${wd}/reads.fa"
zcat -f \
    <(samtools view -h --regions-file "${regions_file}" "${bam}" \
        | samtools fasta -
    ) \
    "${add_fa}" > "${reads}"

bam="${wd}/out.bam"
minimap2 -y -ax lr:hq -I 8G -t 8 "${asm}" "${reads}" \
    | samtools view -u - \
    | samtools sort -o "${bam}"

samtools index "${bam}"

rm -f "${reads}"
