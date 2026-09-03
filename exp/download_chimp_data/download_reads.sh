#/bin/bash

set -euo pipefail

xargs -P 4 -I {} bash -c "wget -nc {}" < chimp_reads.txt
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/033/ERR14368733/ERR14368733.fastq.gz
