#' Preprocessing MS/MS datasets
#'
#' @docType methods
#' @name processMSdata
#' @rdname processMSdata
#'
#' @param proData Proteomics data with rows as gene names and columns as samples
#' @param rnaData Relevant RNAseq datasets for better imputation.
#' @param id The type of gene ids in the input data.
#' @param org Organism.
#' @param k An integer.
#'
#' @return A list including QC figures and intermediate results.
#'
#' @author Wubing Zhang
#'
processMSdata <- function(proData, rnaData = NULL, id = "Symbol", org = "hsa", k = 3){
  message(Sys.time(), " # Transfer ", id, " to Entrez ids ...")
  tmp = TransGeneID(rownames(proData), fromType = id, toType = "Entrez", organism = org)
  idx = is.na(tmp) | duplicated(tmp)
  proData = proData[!idx, ]
  rownames(proData) = tmp[!idx]
  if(!is.null(nrow(rnaData))){
    tmp = TransGeneID(rownames(rnaData), fromType = id, toType = "Entrez", organism = org)
    idx = is.na(tmp) | duplicated(tmp)
    rnaData = rnaData[!idx, ]
    rownames(rnaData) = tmp[!idx]
    rnaData = t(scale(t(rnaData), scale=FALSE))
  }
  ## QC
  message(Sys.time(), " # Quality control ...")
  p1 = countNA(proData)
  gg = data.frame(gene = rownames(proData), NAs = rowSums(is.na(proData)))
  p2 = DensityView(gg[,2,drop=FALSE], main = "NA distribution in gene", xlab = "NA number")
  p2 = p2 + theme(legend.position = "none")
  gg = data.frame(sample = colnames(proData), Detection = colSums(!is.na(proData)))
  p3 = DensityView(gg[,2,drop=FALSE], main = "The number of detected gene", xlab = "Gene")
  p3 = p3 + theme(legend.position = "none")
  ## Normalization
  message(Sys.time(), " # Normalization ...")
  normdata = filterN(proData, 3)
  normdata = normalizeMS(normdata)
  normData = normdata
  if(!is.null(nrow(rnaData))){
    rnaData = as.data.frame(rnaData)
    genes = intersect(rownames(normdata), rownames(rnaData))
    tmp = cbind(normdata[genes, ], rnaData[genes, ])
    rowmax = round(max(rowSums(is.na(tmp))) / ncol(tmp),2)+0.01
    colmax = round(max(colSums(is.na(tmp))) / nrow(tmp),2)+0.01
    imputedata1 = imputeNA(as.matrix(tmp), rowmax = rowmax, colmax = colmax, k = k)
    normdata = rbind(normdata[setdiff(rownames(normdata), genes), ], imputedata1[, 1:ncol(normdata)])
  }
  rowmax = round(max(rowSums(is.na(normdata))) / ncol(normdata),2)+0.01
  colmax = round(max(colSums(is.na(normdata))) / nrow(normdata),2)+0.01
  imputedata = imputeNA(as.matrix(normdata), rowmax = rowmax, colmax = colmax, k = k)
  return(list(p1=p1, p2=p2, p3=p3, rawdata = proData, normData = normData, imputeData = imputedata))
}
