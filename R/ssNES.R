#' single sample gene set enrichment
#' @export
#' @importFrom MAGeCKFlute gsGetter TransGeneID
#' @importFrom GSVA gsva
ssNES <- function(mat, gsType="Complex", method = "gsva"){
  mat = as.matrix(mat)
  gsets = MAGeCKFlute::gsGetter(type = gsType)
  gsets$Symbol = MAGeCKFlute::TransGeneID(gsets$Gene, "Entrez", "Symbol")
  gsets = gsets[!is.na(gsets$Symbol), ]
  gsets_list = unstack(gsets[, c("Symbol", "PathwayName")])
  complexname = gsets[, 2:3]
  complexname = complexname[!duplicated(complexname$PathwayID), ]
  rownames(complexname) = complexname$PathwayID
  gsva_mat <- GSVA::gsva(mat, gsets_list, method = "gsva")
  return(gsva_mat)
}
