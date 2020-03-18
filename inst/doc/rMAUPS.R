## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=FALSE---------------------------------------------------
#  install.packages(c("BiocManager", "devtools"))
#  BiocManager::install("MAGeCKFlute")
#  # OR devtools::install_bitbucket("MAGeCKFlute")
#  devtools::install_bitbucket("liulab/rMAUPS")

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
#  ## Visualize the results on a webpage
#  view("analysis/")

## ----qc--------------------------------------------------------------------
data = normdata[,-1]
p = ViolinView(data, ylab = "Protein abundance")
p = p + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
p
meta = metadata[metadata$Experiment=="experiment1_normdata.csv", -1]
rownames(meta) = meta[,1]
p = pcView(data, meta[colnames(data), 2])
p

## --------------------------------------------------------------------------
data = as.matrix(data)
simulated = data
idx = sample(1:length(simulated), round(0.1*length(simulated)))
simulated[idx] = NA

## ----imputation------------------------------------------------------------
gg = data.frame(gene = rownames(simulated), NAs = rowSums(is.na(simulated)))
p1 = DensityView(gg[,2,drop=FALSE], xlab = "The number of missing value")
p1 + theme(legend.position = "none")
gg = data.frame(sample = colnames(simulated), Detection = colSums(!is.na(simulated)))
p2 = DensityView(gg[,2,drop=FALSE], xlab = "The number of detected gene")
p2 + theme(legend.position = "none")
p3 = BarView(gg, "sample", "Detection", fill = "#8da0cb",
             ylab = "The number of detected gene")
p3 + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))

## --------------------------------------------------------------------------
imputed = imputeNA(simulated)
plot(imputed[idx], data[idx])

## ----dep-------------------------------------------------------------------
## Limma
deres = DEAnalyze(data, meta[,-1], type = "msms", method = "limma")
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

