# .

## Usage
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
