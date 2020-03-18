#' From DEanalyze results to Rabit input
#'
#' @docType methods
#' @name RabitInput
#' @rdname RabitInput
#'
#' @param deres Path to differential expression results or data frame of the DE results.
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
RabitInput <- function(deres, idType = "symbol", org = "hsa"
                       #, filename = NULL
                       ){
  if(!class(deres) %in% c("data.frame", "matrix")){
    diff_gene <- read.table(deres, sep = "\t",stringsAsFactors = FALSE,
                          check.names = FALSE, quote = "", row.names = 1)
  }else{
    diff_gene = deres
  }
  diff_gene <- read.table(deres, sep = "\t", header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE, quote = "", row.names = 1)
  if(tolower(idType) %in% c("symbol", "ensembl")){
    diff_gene$Entrez = TransGeneID(rownames(diff_gene), fromType = idType,
                                   toType = "Entrez", organism = org)
  }else if(tolower(idType) %in% c("uniprot", "refseq")){
    tmp = TransProteinID(rownames(diff_gene), fromType = idType,
                         toType = "symbol", organism = org)
    diff_gene$Entrez = TransGeneID(tmp, fromType = "symbol",
                                   toType = "Entrez", organism = org)
  }else if(tolower(idType) == "entrez"){
    diff_gene$Entrez = rownames(diff_gene)
  }else{
    message("Unrecognized ids: ", paste0(rownames(diff_gene)[1:10],  collapse = ","))
    stop("idType error")
  }

  diff_gene = diff_gene[order(-abs(diff_gene$stat)), ]
  idx = duplicated(diff_gene$Entrez) | is.na(diff_gene$Entrez)
  diff_gene = diff_gene[!idx, ]
  ranklist = diff_gene$stat; names(ranklist) = diff_gene$Entrez
  # if(!is.null(filename))
  #   write.table(ranklist, filename, sep = "\t", quote = FALSE,
  #               row.names = TRUE, col.names = FALSE)
  return(ranklist)
}
