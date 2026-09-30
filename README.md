# Sperm_LR_Amplicon_Recombination_Detection
Detect recombination events in ampliconic regions from LR sperm data.

## Getting started
Clone repo and submodules:
```bash
git clone https://github.com/koisland/sperm_ampl_sv --recursive
cd sperm_ampl_sv
```

Use `pixi` to setup dependencies. Also load `apptainer`.
```bash
module load apptainer
pixi install
```

Compile [`RECallS`](https://github.com/koisland/RECallS) with `rustc` >=1.98.0.
```bash
# To install Rust
# curl --proto '=https' --tlsv1.2 https://sh.rustup.rs -sSf | sh
cargo build --release -m workflow/scripts/RECallS/Cargo.toml
```

## Config
Specify samples for DSA read alignment.
* See https://github.com/logsdon-lab/Snakemake-Aligner for options
```yaml
samples:
  - name: sample_1
    asm_fa: data/asm/sample_1.fa.gz
    # Or uBAM
    reads: [data/reads/sample_1_reads.fastq.gz]
```

If you have a BAM already and a merged diploid assembly:
```yaml
samples:
  - name: sample_1
    asm_fa: data/asm/sample_1.fa.gz
    bam: data/alignments/sample_1.bam
```

### Sections
Currently, this pipeline just does a variety of annotation steps:

#### `call_variants`
* variant calling with `pbsv` and `deepvariant`

#### `liftover`
* assembly to reference alignment with `asm-to-ref-alignment`
* liftover annotation and genes with `ucsc-liftover` and `lifton`

#### `annot`
* repeat masking and segmental duplication annotation with `RepeatMasker` and `biser`
    * TODO: get biser to not eat up memory

#### `call_nahr`
* Call putative intrachromosomal non-allelic homologous recombination events with `RECallS`.

Comment out optional sections (`annot`) you don't want:
```yaml
# annot:
#   ...
```

## Usage
For cluster usage.
```bash
pixi run snakemake \
--configfile config/config.yaml \
--profile workflow/profiles/slurm-executor/ \
-j 50 -np 
```

Otherwise:
```bash
pixi run snakemake \
--configfile config/config.yaml \
-c 50 -np 
```

## Limitations
* If using whole testis, no way to determine germline vs somatic tissues. Previous Schierup study is not reproducible as provided uBAMs do not have methylation information/kinetics.
* Cannot fully validate NAHR events found.
