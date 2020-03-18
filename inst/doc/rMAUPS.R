## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=FALSE---------------------------------------------------
#  install.packages(c("BiocManager", "devtools"))
#  BiocManager::install(c("clusterProfiler", "GSVA", "DESeq2", "limma", "MAGeCKFlute", "msmsTests", "metap", "impute", "ggpubr", "BiocStyle"))
#  devtools::install_github("WubingZhang/rMAUPS")

## ----libs------------------------------------------------------------------
library(MAGeCKFlute)
library(ggplot2)
library(rMAUPS)

## ----preprocess------------------------------------------------------------
## Process multiple datasets in a folder
datapath = system.file("extdata", package = "rMAUPS")
list.files(datapath)
normalizeProteomeDiscoverer(datapath, output = "./", log2 = TRUE)

## Process one dataset
normdata = normalizeProteomeDiscoverer(file.path(datapath, "experiment1_export_proteins.txt"), log2 = TRUE, return = TRUE)
head(normdata)

## ----readData--------------------------------------------------------------
metadata = read.csv(file.path(datapath, "metadata.csv"))
head(metadata)

## ----pipeline, eval=FALSE--------------------------------------------------
#  MAUPSr(system.file("extdata", "metadata.csv", package = "rMAUPS"), outdir = "analysis/")
#  ## Or
#  MAUPSr(metadata, outdir = "analysis/")
#  ## Visualize the results on a shiny app.
#  view("analysis/")
#  
#  # After the shiny app open, please input the path to rMAUPS results, e.g.  "analysis/" here, click `submit`, then all the figure results will be loaded on the webpage. It take seconds to load all the figures, please be patient after clicking `submit`.

## ----simulatedata----------------------------------------------------------
data = as.matrix(normdata[,-1])
meta = metadata[metadata$Experiment=="experiment1_normdata.csv", -1]
rownames(meta) = meta[,1]
simulated = data
idx = sample(1:length(simulated), round(0.1*length(simulated)))
simulated[idx] = NA

## ----qc--------------------------------------------------------------------
qc = ProteomicsQC(simulated, condition = meta[colnames(data), 2], proj.name = "TestQC")
qc$p1
qc$p2
qc$p3
qc$p4
qc$p5
qc$p6
qc$p7

## ----normalize-------------------------------------------------------------
normalized = normalizeProteomics(simulated, norm = "median", log2 = FALSE)

## ----knn-------------------------------------------------------------------
imputed = imputeNA(normalized)
plot(imputed[idx], data[idx])

## ----dep-------------------------------------------------------------------
deres = DEAnalyze(data, meta[,-1], type = "msms", method = "limma")
## Visualize the results using functions from MAGeCKFlute
VolcanoView(deres, "log2FC", "padj", x_cutoff = 0.3, y_cutoff = 0.1)
# Or
deres$logFDR = log10(deres$padj)
ScatterView(deres, x = "log2FC", y = "logFDR", 
            x_cut = c(-0.5,0.5), y_cut = -2, 
            groups = c("bottomleft", "bottomright"), top = 5)

## ----decomplex-------------------------------------------------------------
res = DeComplex(deres)
head(res$deComplex)
res$gobp.p
res$reactome.p
res$gocc.p
res$corum.p

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

