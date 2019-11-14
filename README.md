# rMAUPS (Development)

**rMAUPS** (R project for Model based Analysis of Ubiquitin / Proteasome System) is designed for systematically analyzing proteomic data.

## Installation

```
BiocManager::install(c("clusterProfiler", "GSVA", "DESeq2", "limma", "MAGeCKFlute", "msmsTests", "metap", "impute", "ggpubr", "optparse"))
install.packages("devtools")
devtools::install_github("WubingZhang/rMAUPS")
```

## Availible functions
---
* [x] Normalization.
* [x] Single sample gene set analysis - GSVA, ssGSEA, PC, mean, median, fisher, Edgington, ...
* [x] Differential expression analysis - limma, DEqMS, DESeq2, ... 
* [x] Pratt and MEME API.
* [x] Protein ID - Gene ID conversion.
* [] Imputation. (Incorporate functions in imputeLCMD)
* [] Incoporate rMAUPS into MAUPS snakemake pipeline.

---
