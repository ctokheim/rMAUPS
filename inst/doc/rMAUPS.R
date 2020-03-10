## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=FALSE-------------------------------------------------
#  install.packages(c("BiocManager", "devtools", "ggplot2"))
#  BiocManager::install("MAGeCKFlute")
#  devtools::install_bitbucket("liulab/rMAUPS")

## ----libs----------------------------------------------------------------
library(rMAUPS)
library(MAGeCKFlute)
library(ggplot2)

## ----preprocess, eval=FALSE----------------------------------------------
#  folder = "examples/rawdata/"
#  normalizeProteomeDiscoverer(folder, log2 = TRUE)

## ----pipeline, eval=FALSE------------------------------------------------
#  ## Download test data from bitbucket
#  setwd("testdata")
#  MAUPSr(metadata = "metadata.csv", outdir = "../analysis/")

## ----readData------------------------------------------------------------
data("testdata")
head(testdata)
meta = data.frame(row.names = colnames(testdata),
                  Condition = gsub(".*\\.", "", colnames(testdata)),
                  stringsAsFactors = FALSE)
meta

## ----qc------------------------------------------------------------------
p = ViolinView(testdata, ylab = "Protein abundance")
p = p + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
p

p = pcView(testdata, meta$Condition)
p

## ----imputation----------------------------------------------------------
## Randomly assign 10% values to be NA
testdata = as.matrix(testdata)
simulated = testdata
idx = sample(1:length(simulated), round(0.1*length(simulated)))
simulated[idx] = NA
## Impute missing values using KNN
imputed = imputeNA(simulated)
plot(imputed[idx], testdata[idx])

## ----dep-----------------------------------------------------------------
## Limma
deres = DEAnalyze(testdata, meta, type = "msms", method = "limma")
VolcanoView(deres, "log2FC", "padj", x_cutoff = 0.3, y_cutoff = 0.1)

## ----decomplex-----------------------------------------------------------
res = deComplex(deres)
head(res$deComplex)
res$gobp.p
res$reactome.p
res$gocc.p
res$corum.p

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

