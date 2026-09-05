#/bin/bash

xargs -P 4 -I {} bash -c "wget {}" < chimp_asm.txt
# remove api stuff from filename
# remove comments and '|' from assembly sequence names
rename "?download=true" ".fa" GCA*
find *.fa | \
xargs -P 4 -I {} \
    bash -c "seqkit replace {} -p "\s.+" | \
    sed 's/|/_/g' | \
    bgzip > {}.gz && \
    samtools faidx {}.gz"
