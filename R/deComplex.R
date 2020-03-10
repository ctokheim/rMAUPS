##' Convert gene level differences to pathway level or complex level differences.
##'
##' @param depres normalized count matrix; rows are all genes in the signature that
##' shall be summarized into one score; columns are samples
##' @param lfc Gene sets.
##' @param pval ("PC", default), Pearson, ssGSEA or mean (other value). fisher, stouffer
##' @return A list.
##' @author Wubing Zhang
##' @export
##'
deComplex <- function(depres, lfc = "log2FC", pval = "pvalue"){
  ## GOCC, CORUM, Pathway
  gsets = gsGetter(type = "GOBP+GOCC+CORUM+REACTOME", limit = c(0,200))
  gsets$Symbol = TransGeneID(gsets$Gene, "Entrez", "Symbol")
  gsets = gsets[gsets$Symbol %in% rownames(depres), ]
  genesets = unstack(gsets[,c(4,2)])

  fisherp = gsScore(depres[,pval, drop=FALSE], gset = genesets, fun = "fisher")
  stoufferZ = gsScore(depres[,lfc, drop=FALSE], gset = genesets, fun = "stouffer")
  merged_deres = cbind(Zscore = stoufferZ, pvalue = fisherp)
  merged_deres = as.data.frame(merged_deres, stringsAsFactors = FALSE)
  colnames(merged_deres) = c("Zscore", "pvalue")
  tmp = gsets[!duplicated(gsets$PathwayID), ]
  rownames(tmp) = tmp$PathwayID
  merged_deres$Description = tmp[rownames(merged_deres), 3]
  merged_deres = merged_deres[, c(3,1,2)]

  ## Visualize differential expressed complexes and pathways
  merged_deres$logP = -log10(merged_deres$pvalue)
  tmp = merged_deres[grepl("GOBP", rownames(merged_deres)), ]
  p1 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 5, main = "GOBP", ylab = "-log10(p.value)")
  p1 = p1 + theme(legend.position = "none")
  tmp = merged_deres[grepl("REACTOME", rownames(merged_deres)), ]
  p2 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 5, main = "REACTOME", ylab = "-log10(p.value)")
  p2 = p2 + theme(legend.position = "none")
  tmp = merged_deres[grepl("GOCC", rownames(merged_deres)), ]
  p3 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 5, main = "GOCC", ylab = "-log10(p.value)")
  p3 = p3 + theme(legend.position = "none")
  tmp = merged_deres[grepl("CORUM", rownames(merged_deres)), ]
  p4 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 5, main = "CORUM", ylab = "-log10(p.value)")
  p4 = p4 + theme(legend.position = "none")
  return(list(deComplex = merged_deres, gobp.p = p1, reactome.p = p2,
              gocc.p = p3, corum.p = p4))
}
