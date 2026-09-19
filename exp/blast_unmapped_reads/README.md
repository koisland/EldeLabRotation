# Purpose
Unclear whether or not recombination events are merely unmapped reads.

Test via:
1. Pull unmapped reads from read alignment
1. Subset reference assembly to region of interest
1. Repeatmask region and generate masked region of interest
1. Run blast to see if any homology

```bash
pushd exp/blast_unmapped_reads
bash run.sh input.bam input_ref.fa chr7:1-1000
```

## Outcome
Nothing useful. Few reads (100s) and no large hits for gene of interest.
