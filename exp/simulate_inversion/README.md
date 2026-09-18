# Purpose
Generate simulated reads to observe read alignment patterns in ampliconic regions via:
* Existing haplotypes
* Simulated inversions/duplications from donor haplotype

## Existing haplotypes
Simulate PacBio HiFi reads of existing haplotypes.
* using `pbsim`.
```bash
bash exp/simulate_inversion/existing_haps_generate_sim.sh \
    exp/simulate_inversion/representative_tcaf.fa.gz \
    exp/simulate_inversion/
```

Will create directories:
```
exp/simulate_inversion/fasta/
└── split_NA21093#2#CM089609.1__144539366-145195716.fa
exp/simulate_inversion/sim/
├── ...
├── NA21093#2#CM089609.1:144539366-145195716_0001.fq.gz
├── NA21093#2#CM089609.1:144539366-145195716_0001.maf.gz
└── NA21093#2#CM089609.1:144539366-145195716_0001.ref
```

Align an existing haplotype as simulated hifi reads to a region of interest.
* Using `minimap2` and the `map-hifi` preset
```bash
sample=""
bash exp/simulate_inversion/existing_haps_align_rgn.sh \
    "results/align/${sample}.bam" \
    "results/align/${sample}.fa" \
    "exp/simulate_inversion/${sample}_TCAF_regions.bed" \
    <(zcat NA21093#2#CM089609.1:144539366-145195716_0001.fq.gz)
```

## Simulated inversions/duplications from donor haplotypes
Drawing dotplots with `python` and `minimap2`
* We use recommended parameters (https://github.com/lh3/minimap2/issues/106) for dotplot visualization.

```bash
python exp/simulate_inversion/denovo_inv_create.py \
-f "exp/simulate_inversion/fasta/split_HG00272#2#CM094217.1__145896774-146653873.fa" \
-s 42
```
