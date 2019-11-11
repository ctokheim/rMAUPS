#' single sample gene set enrichment
#' @export
ssNES <- function(mat, gsType="Complex", method = "gsva"){
  require(MAGeCKFlute)
  require(GSVA)
  mat = as.matrix(mat)
  gsets = gsGetter(type = gsType)
  gsets$Symbol = TransGeneID(gsets$Gene, "Entrez", "Symbol")
  gsets = gsets[!is.na(gsets$Symbol), ]
  gsets_list = unstack(gsets[, c("Symbol", "PathwayName")])
  complexname = gsets[, 2:3]
  complexname = complexname[!duplicated(complexname$PathwayID), ]
  rownames(complexname) = complexname$PathwayID
  gsva_mat <- gsva(mat, gsets_list, method = "gsva")
  return(gsva_mat)
}
