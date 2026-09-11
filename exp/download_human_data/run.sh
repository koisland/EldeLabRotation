#!/bin/bash

set -euo pipefail

outdir="data/annot/human"
mkdir -p "${outdir}"

# assembly
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -P "${outdir}"
# genes
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz -P "${outdir}"
gunzip "${outdir}/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz" 
# segdups
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/CHM13/bbi/segDups_2024.bb -P "${outdir}"
bigbedtobed "${outdir}/segDups_2024.bb" "${outdir}/segDups_2024.bed"
