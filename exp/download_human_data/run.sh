#!/bin/bash

set -euo pipefail

outdir="data/annot/human"
mkdir -p "${outdir}"

# assembly
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -P "${outdir}"

# genes
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/browser/bbi/chm13v2.0_RefSeq_Liftoff_v5.2.bb -P "${outdir}"
bigbedtobed "${outdir}/segDups_2024.bb" "${outdir}/chm13v2.0_RefSeq_Liftoff_v5.2.bed"
awk -v OFS="\t" '{ print $1, $2, $3, $4","$19, $5, $6, $7, $8, "0,0,0", $10, $11, $12 }' "${outdir}/chm13v2.0_RefSeq_Liftoff_v5.2.bed" > "${outdir}/chm13v2.0_RefSeq_Liftoff_v5.2_bed12.bed"

# segdups
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/CHM13/bbi/segDups_2024.bb -P "${outdir}"
bigbedtobed "${outdir}/segDups_2024.bb" "${outdir}/segDups_2024.bed"
