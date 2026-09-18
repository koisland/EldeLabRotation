# Purpose
Generate simulated reads to observe read alignment patterns in ampliconic regions via:
* Existing haplotypes
* Simulated inversions/duplications from donor haplotype

## Existing haplotypes
Split existing haplotypes into non-overlapping intervals (20 kbp reads) and add IT:Z tag to header.
```bash
bash exp/simulate_inversion/exisitng_haps_generate_sim.sh \
    exp/simulate_inversion/representative_tcaf.fa.gz \
    exp/simulate_inversion/
```

Align existing haplotypes as 20 kbp reads to region of interest.
```bash
bash exp/simulate_inversion/existing_haps_align_rgn.sh \
    results/align/0068.bam \
    results/align/0068.fa \
    exp/simulate_inversion/0068_TCAF_regions.bed \
    <(cat exp/simulate_inversion/*.fa)
```

```bash
bash exp/simulate_inversion/existing_haps_align_rgn.sh \
    results/align/0068.bam \
    results/align/0068.fa \
    exp/simulate_inversion/0068_TCAF_regions.bed \
    <(zcat exp/simulate_inversion/sim/HG03041#2#CM088732.1:144576079-144943903_0001.fq.gz)
```

## Simulated inversions/duplications from donor haplotypes
Drawing dotplots with tenten and minimap2
* https://github.com/ocxtal/tenten
```bash
$ cargo install --git https://github.com/ocxtal/tenten.git
tenten -s exp/simulate_inversion/HG00272#2#CM094217.1:145896774-146653873.fa -b 500
```

```bash
```