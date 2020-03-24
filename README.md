# rMAUPS

Advances in proteomic profiling have enabled the study of protein regulation and degradation in diseases such as cancer. However, there are major computational challenges, such as how to perform quality control (QC) and normalization, how to identify differential protein abundance at multiple-levels of resolution (e.g. protein complexes), and how to integrate data with other omic technologies. Here, we developed a computational analysis pipeline Model-based Analysis of the Ubiquitin-Proteasome System using R (**rMAUPS**), which performs computational analysis of proteomics data efficiently and effectively. rMAUPS comprises four significant modules, including quality control (QC), imputation, differential analysis, and integrative analysis.

## Installation
Installing the package in a fresh R environment may take a long time. It may fail because of some issues. You can check the possible issues and solutions from https://github.com/WubingZhang/rMAUPS/issues/3, or post a new issue there.

## Installation using R
One of the dependencies require `libnetcdf`, you can install it using `conda install -c anaconda libnetcdf` in the bash command line. 
```
$ conda install -c anaconda libnetcdf
```
Then you can install all the required R packages in R.
```
> install.packages(c("devtools", "BiocManager"), repos = "https://cloud.r-project.org")
> BiocManager::install(c("DESeq2", "limma", "impute", "biomaRt", "GO.db", "msmsTests"))
> devtools::install_github("WubingZhang/rMAUPS")
```

## Installation using bash command line
```
$ conda install -c anaconda libnetcdf
$ conda install -c conda-forge r-devtools r-biocmanager r-ggpubr r-metap
$ conda install -c bioconda bioconductor-gsva bioconductor-biocstyle
$ Rscript -e 'BiocManager::install(c("DESeq2", "limma"), quiet=T)'
$ Rscript -e 'BiocManager::install(c("msmsTests"))'
$ Rscript -e 'BiocManager::install(c("impute"), quiet=T)'
$ Rscript -e 'BiocManager::install(c("biomaRt", "GO.db"), quiet=T)'
$ Rscript -e 'devtools::install_github("WubingZhang/rMAUPS")'
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

