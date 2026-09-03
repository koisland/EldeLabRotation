# Data
Get chimp reads from paper: [*Long-read sequencing of primate testis and human sperm allows identification of recombination events in individuals*](https://www.nature.com/articles/s41467-025-65248-3)

Start with chimpanzee as example case since human data is pending.

All data under ENA accession [PRJEB77177](https://www.ebi.ac.uk/ena/browser/view/PRJEB77177)

## Assemblies
See example: https://www.ebi.ac.uk/ena/browser/view/CBCUDL010000000?show=sample-attributes

* CT15: GCA_965153145.1
* CT22: GCA_965153165.1
* CT28: GCA_965153115.1
* CT32: GCA_965153035.1
    * Not explicitly shown sample name in sample attributes but can infer since last chimp

## Reads
See https://www.ebi.ac.uk/ena/browser/view/PRJEB77177?show=reads

* CT15: ERR14376286
* CT22: ERR14376368
* CT28: ERR14376446
* CT32: ERR14368733

## Run
```bash
bash download_asm.sh
# remove api stuff from filename
rename "?download=true" ".fa.gz" GCA*

bash download_reads.sh
```
