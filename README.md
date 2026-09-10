# Sperm_LR_Amplicon_Recombination_Detection
Detect recombination events in ampliconic regions from LR sperm data.

## Usage
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

Then run:
```bash
pixi run snakemake \
--configfile config/config.yaml \
--profile workflow/profiles/slurm-executor/ \
-j 50 -np 
```

## Limitations
* If using whole testis, no way to determine germline vs somatic tissues. Previous Schierup study is not reproducible as provided uBAMs do not have methylation information/kinetics.
