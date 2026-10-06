#!/bin/bash

set -eu

wd=$(dirname $0)

python workflow/scripts/visualize_events.py \
-i <(
    cat \
    <(awk -v OFS="\t" '{ print $1, $2, $3, "Inversion signal", 0, ".", $2, $3, "0,0,255" }' results/call_nahr/0068/calls_inv.bed) \
    <(awk -v OFS="\t" '{ print $1, $2, $3, "Deletion signal", 0, ".", $2, $3, "255,0,0" }' results/call_nahr/0068/calls_del.bed)
) \
-a <(zcat results/liftover/CHM13/annot/0068_sd.bed.gz | cut -f1-9) \
-f results/align/0068.fa.fai \
-c <(awk '$3-$2 > 100000' results/annot/longdust/0068_lcr.bed)  \
-o "${wd}/ideogram"
