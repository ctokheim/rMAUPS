# rMAUPS (Development)

**rMAUPS** (R project for Model based Analysis of Ubiquitin / Proteasome System) is designed for systematically analyzing proteomic data.

## Installation

```
BiocManager::install(c("clusterProfiler", "GSVA", "DESeq2", "limma", "MAGeCKFlute", "msmsTests", "metap", "impute", "ggpubr", "BiocStyle"))
install.packages("devtools")
devtools::install_github("WubingZhang/rMAUPS")
```

## Tutorial
For detail documentation, please visit http://cistrome.org/~wubing/rMAUPS.html or run the following command line.

```
library(rMAUPS)
vignette("rMAUPS")
```


## Contacts

* Wubing Zhang (watson5bzhang@gmail.com)
* Collin Tokheim (ctokheim@ds.dfci.harvard.edu)

