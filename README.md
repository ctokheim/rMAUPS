# rMAUPS

Advances in proteomic profiling have enabled the study of protein regulation and degradation in diseases such as cancer. However, there are major computational challenges, such as how to perform quality control (QC) and normalization, how to identify differential protein abundance at multiple-levels of resolution (e.g. protein complexes), and how to integrate data with other omic technologies. Here, we developed a computational analysis pipeline Model-based Analysis of the Ubiquitin-Proteasome System using R (**rMAUPS**), which performs computational analysis of proteomics data efficiently and effectively. rMAUPS comprises four significant modules, including quality control (QC), imputation, differential analysis, and integrative analysis.

## Installation
Installing the package in a fresh R environment may take a long time. It may fail because of some issues. You can check the possible issues and solutions from https://github.com/WubingZhang/rMAUPS/issues/3, or post a new issue there.

### Install rMAUPS using bash command line (recommended)
```
conda install -c bioconda bioconductor-gsva bioconductor-deseq2 bioconductor-limma bioconductor-msmstests bioconductor-impute bioconductor-biomart bioconductor-biocstyle
Rscript -e 'devtools::install_github("WubingZhang/rMAUPS")'
```
### Install rMAUPS using R command line
```
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("GSVA", "DESeq2", "limma", "msmsTests", "metap", "impute", "BiocStyle", "biomaRt")) 
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

