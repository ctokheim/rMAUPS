#' single sample gene set enrichment
#'
#' @docType methods
#' @name DeComplex
#' @rdname DeComplex
#'
#' @param mat A matrix.
#' @param gsType Molecular signatures for testing, available datasets include Pathway
#' (KEGG, REACTOME, C2_CP), GO (GOBP, GOCC, GOMF), MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA),
#' C2_CGP), C3 (C3_MIR, C3_TFT), C4, C6, C7, HALLMARK) and Complex (CORUM). Any combination of them
#' are also accessible (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param method Method to employ in the estimation of gene-set enrichment scores per sample.
#' By default this is set to gsva (Hänzelmann et al, 2013) and other options are ssgsea
#' (Barbie et al, 2009), zscore (Lee et al, 2008) or plage (Tomfohr et al, 2005). The latter two
#' standardize first expression profiles into z-scores over the samples and, in the case of zscore,
#' it combines them together as their sum divided by the square-root of the size of the gene set,
#' while in the case of plage they are used to calculate the singular value decomposition (SVD) over
#' the genes in the gene set and use the coefficients of the first right-singular vector as pathway
#' activity profile.
#'
#' @importFrom MAGeCKFlute gsGetter TransGeneID
#' @importFrom GSVA gsva
#' @export
ssNES <- function(mat, gsType="Complex", method = "gsva"){
  mat = as.matrix(mat)
  gsets = MAGeCKFlute::gsGetter(type = gsType)
  gsets$Symbol = MAGeCKFlute::TransGeneID(gsets$Gene, "Entrez", "Symbol")
  gsets = gsets[!is.na(gsets$Symbol), ]
  gsets_list = unstack(gsets[, c("Symbol", "PathwayName")])
  complexname = gsets[, 2:3]
  complexname = complexname[!duplicated(complexname$PathwayID), ]
  rownames(complexname) = complexname$PathwayID
  gsva_mat <- GSVA::gsva(mat, gsets_list, method = method)
  return(gsva_mat)
}
