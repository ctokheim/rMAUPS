#' From DEanalyze results to Rabit input
#'
#' @docType methods
#' @name RabitInput
#' @rdname RabitInput
#'
#' @param DEpath Path to differential expression results.
#' @param idType One of "symbol", "ensembl", "uniprot", and "refseq".
#' the rownames should match colnames in obj, and the first column should be Condition.
#' @param type "Array", "RNASeq" or "msms", only needed when obj is matrix like object.
#' @param method Differential expression analysis method, e.g. limma, DESeq2, GFOLD,
#' glm.pois, glm.qlll, and glm.nb.
#' @param paired Boolean, specifying whether perform paired comparison.
#' @param app.dir The path to application (e.g. GFOLD).
#'
#' @return An ExpressionSet instance.
#' @seealso \code{\link{ExpressionSet-class}}
#'
#' @author Wubing Zhang
#'
#' @import limma DESeq2 msmsTests Biobase
#' @export
RabitInput <- function(DEpath, idType = "symbol", org = "hsa", outdir = "./"){
  diff_gene <- read.table(DEpath, sep = "\t", header = TRUE, check.names = FALSE,
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
  write.table(ranklist, paste0(outdir, "/RabitInput_", DEpath),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
}
