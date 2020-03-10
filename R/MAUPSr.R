##' rMAUPS pipeline - QC, differential analysis, integrative analysis.
##'
##' @param metadata File path or data frame of the meta data, including columns of
##' Experiment, Sample, Condition, and multiple comparisons.
##' @param qc Boolean, specifying whether perform the quanlity control.
##' @param outdir Output directory.
##' @param type Could be "msms", "RNASeq", or "Arrary".
##'
##' @author Wubing Zhang
##' @return No return value. Output multiple files into local folder.
##' @export
##'
MAUPSr <- function(metadata, qc = TRUE, outdir = "./", type = "msms"){
  options(stringsAsFactors = FALSE)
  message(format(Sys.time(), "%b-%d-%Y %X "), "Reading the metadata ...")
  if(length(metadata)==1 && file.exists(metadata)){# Read meta file
    if(grepl(".csv$", metadata))
      metadata = read.csv(metadata, header = TRUE)
    else
      metadata = read.table(metadata, sep = "\t", header = TRUE)
  }

  message(format(Sys.time(), "%b-%d-%Y %X "), "Analyzing experiments one by one ...")
  ## Analyze each data separately
  for(experiment in unique(metadata[,1])){
    message("\t", format(Sys.time(), "%b-%d-%Y %X "), experiment, " ...")
    data = read.csv(experiment, header = TRUE, row.names = 1)
    meta = metadata[metadata[,1]==experiment, ]
    rownames(meta) = meta[,2]
    meta = meta[, colSums(is.na(meta))<nrow(meta)]
    proj.name = gsub(".*\\/|_normdata.*|\\..*", "", experiment)
    if(qc) plist = qcProteomics(data, condition = meta[colnames(data), 3],
                                proj.name = proj.name, outdir = outdir)
    if(sum(is.na(data))>0){
      rowmax = round(max(rowSums(is.na(data))) / ncol(data),2)+0.01
      colmax = round(max(colSums(is.na(data))) / nrow(data),2)+0.01
      data = filterN(data)
      imputedata = imputeNA(as.matrix(data), rowmax = rowmax, colmax = colmax, k = 5)
      tmpfile = paste0(outdir,"/",basename(gsub("\\..*", "_imputed.csv", experiment)))
      write.csv(imputedata, tmpfile, row.names = TRUE, quote = FALSE)
      data = imputeddata
    }
    comparisons = grep("comparison", colnames(meta),
                       ignore.case = TRUE, value = TRUE)
    for(comp in comparisons){
      prefix = paste0(proj.name, ".", comp)
      SA = meta[!is.na(meta[,comp]), comp, drop = FALSE]
      message("\t", format(Sys.time(), "%b-%d-%Y %X "), comp, " ...")
      method = ifelse(grepl("RNA", type), "DESeq2", "limma")
      deres_p = DEAnalyze(data, SA, type = type, method = method)
      write.csv(deres_p, paste0(outdir,"/",prefix, "_dep.csv"),
                row.names = TRUE, quote = FALSE)
      deres_p$logP = -log(deres_p$pvalue)
      p1 = ScatterView(deres_p, x = "log2FC", y = "logP",
                       model = "volcano", x_cut = c(-0.2,0.2), force = 5,
                       top = 5, ylab = "-log10(p.value)")
      p1 = p1 + theme(legend.position = "none")
      ggsave(paste0(outdir,"/",prefix, "_dep_volcano.png"),
             p1, width = 6, height = 5)

      res = deComplex(deres_p)
      write.csv(res$deComplex, paste0(outdir,"/",prefix, "_dePathway.csv"),
                row.names = TRUE, quote = FALSE)
      ggsave(paste0(outdir,"/",prefix, "_deBP_volcano.png"),
             res$gobp.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/",prefix, "_deREACTOME_volcano.png"),
             res$reactome.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/",prefix, "_deCC_volcano.png"),
             res$gocc.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/",prefix, "_deCC_volcano.png"),
             res$corum.p, width = 6, height = 5)
    }
  }
  ## Merge the same comparisons in different experiments
  message(format(Sys.time(), "%b-%d-%Y %X "),
          "Merge the same comparisons in different experiments ...")
  comparisons = grep("comparison", colnames(metadata), ignore.case = TRUE, value = TRUE)
  for(comp in comparisons){
    experiments = unique(metadata[!is.na(metadata[,comp]),1])
    if(length(experiments)>1){
      message("\t", format(Sys.time(), "%b-%d-%Y %X "), comp, " ...")
      proj.names = gsub(".*\\/|_normdata.*|\\..*", "", experiments)
      prefix = paste0(proj.names, ".", comp)
      DEPs = paste0(outdir,"/",prefix, "_dep.csv")
      summary = data.frame()
      for(r in DEPs){
        tmp = read.csv(r, row.names = 1, header = TRUE)
        proteins = unique(c(rownames(summary), rownames(tmp)))
        summary = cbind(summary[proteins,], tmp[proteins, c(1,4)])
        rownames(summary) = proteins
      }
      colnames(summary) = paste0(rep(proj.names,each=2), ".", rep(c("log2FC", "pvalue"),length(DEPs)))
      mergedDep = t(apply(summary, 1, function(x){
        lfc = x[seq(1, length(x),2)]; pval = x[seq(2, length(x),2)]
        lfc = lfc[!is.na(lfc)]; pval = pval[!is.na(pval)]
        if(length(lfc)==0) return(c(NA, NA))
        if(length(lfc)==1) return(c(lfc, pval))
        if(length(lfc)>1) c(sum(lfc)/sqrt(length(lfc)), metap::sumlog(pval)$p)
      }))
      mergedDep = as.data.frame(mergedDep, stringsAsFactors = FALSE)
      colnames(mergedDep) = c("Zscore", "pvalue")
      write.csv(mergedDep, paste0(outdir,"/", comp, "_merged.dep.csv"), quote = FALSE)

      res = deComplex(mergedDep, lfc = "Zscore")
      write.csv(res$deComplex, paste0(outdir,"/", comp, "_merged.dePathway.csv"),
                row.names = TRUE, quote = FALSE)
      ggsave(paste0(outdir,"/", comp, "_merged_deBP_volcano.png"),
             res$gobp.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/", comp, "_merged_deREACTOME_volcano.png"),
             res$reactome.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/", comp, "_merged_deCC_volcano.png"),
             res$gocc.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/", comp, "_merged_deCORUM_volcano.png"),
             res$corum.p, width = 6, height = 5)
    }
  }
}
