# rMAUPS (Development)

**rMAUPS** (R project for Model based Analysis of Ubiquitin / Proteasome System) is designed for systematically analyzing proteomic data.

## Installation

```
BiocManager::install(c("clusterProfiler", "GSVA", "DESeq2", "limma", "MAGeCKFlute", "msmsTests", "metap", "impute", "ggpubr", "optparse"))
install.packages("devtools")
devtools::install_github("WubingZhang/rMAUPS")
```

## Visualize tutorial
```
vignette("rMAUPS")
```

## Availible functions
---
* [x] Data quality control and normalization.
* [x] Differential abundance analysis - limma, DEqMS, DESeq2, ... 
* [x] Protein complex level analysis
* [x] Pratt and MEME API.
* [x] Protein ID - Gene ID conversion.
* [x] Imputation. (Todo: Incorporate functions in imputeLCMD)
* [] Incoporate rMAUPS into MAUPS snakemake pipeline.
---
