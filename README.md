# rMAUPS

Advances in proteomic profiling have enabled the study of protein regulation and degradation in diseases such as cancer. However, there are major computational challenges, such as how to perform quality control (QC) and normalization, how to identify differential protein abundance at multiple-levels of resolution (e.g. protein complexes), and how to integrate data with other omic technologies. Here, we developed a computational analysis pipeline Model-based Analysis of the Ubiquitin-Proteasome System using R (**rMAUPS**), which performs computational analysis of proteomics data efficiently and effectively. rMAUPS comprises four significant modules, including quality control (QC), imputation, differential analysis, and integrative analysis.

## Installation
Installing the package in a fresh environment may take a long time. It may fail because of some issues. You can check the possible issues and solutions from https://github.com/WubingZhang/rMAUPS/issues/3, or post a new issue there.

### Prerequisites
To install rMAUPS, you have to first install conda following the document (https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html#install-macos-silent).  The fast way is shown as follows:  

```
# Installing conda on macOS
wget wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Installing conda on Linux
wget wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Installing conda on Windows
# 1. Download the installer (.exe file) from https://conda.io/miniconda.html
# 2. Install miniconda by double clicking the file
# 3. Test your installation. In your terminal window or Anaconda Prompt, run the command conda list.
```

You should also install r, r-recommended and libnetcdf (required by rMAUPS) in your conda environment.
```
$ conda install -c r r r-recommended r-markdown
$ conda install -c anaconda libnetcdf
```

### Installation of rMAUPS using R
rMAUPS is a R package released from github, so you can install the it using R functions.
```
> install.packages(c("devtools", "BiocManager"), repos = "https://cloud.r-project.org")
# Install dependencies
> BiocManager::install(c("GSVA", "DESeq2", "limma", "impute", "biomaRt", "GO.db", "msmsTests", "BiocStyle"))
# Install rMAUPS from github
> devtools::install_github("WubingZhang/rMAUPS")
```

### Installation using bash command line
```
# Install dependencies using conda
$ conda install -c conda-forge r-devtools r-biocmanager r-ggpubr r-metap
# Installation of some packages from Bioconductor using conda always generates conflicts, so please install these dependencies using Rscript.
$ Rscript -e 'BiocManager::install(c("GSVA", "DESeq2", "limma", "impute", "biomaRt", "GO.db", "msmsTests", "BiocStyle"))'
$ Rscript -e 'devtools::install_github("WubingZhang/rMAUPS")'
```

## Tutorial
For detail documentation, please visit http://cistrome.org/~wubing/rMAUPS.html or run the following command line.

```
> library(rMAUPS)
> vignette("rMAUPS")
```


## Contacts

* Wubing Zhang (watson5bzhang@gmail.com)
* Collin Tokheim (ctokheim@ds.dfci.harvard.edu)

