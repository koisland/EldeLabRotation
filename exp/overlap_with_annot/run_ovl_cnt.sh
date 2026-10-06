#!/bin/bash

set -eu

wd=$(dirname $0)

for ftype in "del" "inv"; do
    bedtools intersect \
    -a <(sort -k1,1 -k2,2n -k3,3 -k4,4 -k5,5n -u /scratch/ucgd/lustre-labs/vollger/users/Keith/EldeLabRotation/results/call_nahr/0068/calls_${ftype}.bed) \
    -b <(zcat /scratch/ucgd/lustre-labs/vollger/users/Keith/EldeLabRotation/results/liftover/CHM13/annot/0068_sd.bed.gz | cut -f1-3 | grep chr) -loj \
    > "${wd}/ovl_${ftype}_sd.bed"

    bedtools intersect \
    -a /scratch/ucgd/lustre-labs/vollger/users/Keith/EldeLabRotation/results/call_nahr/0068/calls_${ftype}.bed \
    -b <(zcat /scratch/ucgd/lustre-labs/vollger/users/Keith/EldeLabRotation/results/liftover/CHM13/annot/0068_genes.bed.gz | cut -f1-4 | grep chr) -loj \
    > "${wd}/ovl_${ftype}_genes.bed"
done
