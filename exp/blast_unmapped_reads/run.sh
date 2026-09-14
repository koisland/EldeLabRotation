#!/bin/bash

set -euo pipefail

# Update pixi.toml
# repeatmasker = "==4.1.2.p1"
# blast = { version = ">=2.16.0,<3", channel = "bioconda" }

wd=$(dirname $0)
bam="${1}"
ref="${2}"
rgn="${3}"

unmap_bam="${wd}/unmapped.bam"
unmap_fa="${wd}/unmapped.fa"

# Get unmapped reads
samtools view -f 4 -@ 16 -o "/${unmap_bam}" "${bam}"
samtools fasta "${unmap_bam}" > "${unmap_fa}"

# Generate library
# * run repeatmasker and softmask regions
# * then convert to blast mask
rgn_subset="${wd}/${rgn}.fa"
samtools faidx "${ref}" "${rgn}" > "${rgn_subset}"

rm_outdir="${wd}/${rgn}_rm"
RepeatMasker -xsmall -species human -pa 4 -engine rmblast -dir "${rm_outdir}" "${rgn_subset}"

masked_rgn_fa="${rm_outdir}/${rgn}.fa.masked"
out_blast_mask="${wd}/${rgn}.fa.asnb"
out_blast_mask_db="${wd}/${rgn}.fa.asnb.db"
convert2blastmask \
    -in "${masked_rgn_fa}" \
    -parse_seqids \
    -masking_algorithm repeat \
    -masking_options "repeatmasker, default" \
    -outfmt maskinfo_asn1_bin \
    -out "${out_blast_mask}"

makeblastdb \
    -in "${masked_rgn_fa}" \
    -dbtype nucl \
    -parse_seqids \
    -mask_data "${out_blast_mask}" \
    -out "${out_blast_mask_db}" \
    -title "${rgn}"

# query_id	subject_id	per_identity	aln_length	mismatches	gap_openings	q_start	q_end	s_start	s_end	e-value	bit_score
# https://rnnh.github.io/bioinfo-notebook/docs/blast.html#blast--outfmt-6-results
blastn -query "${unmap_fa}" -db ${out_blast_mask_db} -out "${wd}/${rgn}_blast.out" -outfmt 6 -evalue 1e-30
