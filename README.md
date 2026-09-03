# Sperm_LR_Amplicon_Recombination_Detection
Detect recombination events in ampliconic regions from LR sperm data.

## Usage
Clone repo and submodules:
```bash
git clone https://github.com/koisland/sperm_ampl_sv --recursive
cd sperm_ampl_sv
```

Use `pixi` to setup dependencies.
```bash
pixi install
```

Then run:
```bash
pixi run snakemake \
--configfile config/config.yaml \
--profile workflow/profiles/slurm-executor/ \
-j 50 -np 
```
