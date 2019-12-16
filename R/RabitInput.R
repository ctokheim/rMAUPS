#' From DEanalyze results to Rabit input
#'
#' @docType methods
#' @name RabitInput
#' @rdname RabitInput
#'
#' @param DExprSet Data frame of the DE results.
#' @param idType One of "symbol", "ensembl", "uniprot", and "refseq".
#' the rownames should match colnames in obj, and the first column should be Condition.
#' @param org hsa or mmu.
#'
#' @return A weighted gene list.
#'
#' @author Wubing Zhang
#'
#' @import MAGeCKFlute
#' @export
RabitInput <- function(DExprSet, idType = "symbol", org = "hsa"){
  if(tolower(idType) %in% c("symbol", "ensembl")){
    DExprSet$Entrez = TransGeneID(rownames(DExprSet), fromType = idType,
                                  toType = "Entrez", organism = org)
  }else if(tolower(idType) %in% c("uniprot", "refseq")){
    tmp = TransProteinID(rownames(DExprSet), fromType = idType,
                         toType = "symbol", organism = org)
    DExprSet$Entrez = TransGeneID(tmp, fromType = "symbol",
                                  toType = "Entrez", organism = org)
  }else if(tolower(idType) == "entrez"){
    DExprSet$Entrez = rownames(DExprSet)
  }else{
    message("Unrecognized ids: ", paste0(rownames(DExprSet)[1:10],  collapse = ","))
    stop("idType error")
  }

  DExprSet = DExprSet[order(-abs(DExprSet$stat)), ]
  idx = duplicated(DExprSet$Entrez) | is.na(DExprSet$Entrez)
  DExprSet = DExprSet[!idx, ]
  ranklist = DExprSet$stat; names(ranklist) = DExprSet$Entrez
  return(ranklist)
}
