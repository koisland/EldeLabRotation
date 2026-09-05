
#!/bin/bash

set -euo pipefail

mkdir -p annot

# segdups
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0.SD_v1.0.bb -P annot
# gene
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/gene/mPanTro3_refGene.bb -P annot
