#' PreprocessCPTAC
#' @docType methods
#' @name PreprocessCPTAC
#' @rdname PreprocessCPTAC
#' @param expr Protein abundance data or Phosphorylation site abundance data.
#' @param selection Unshared.Log.Ratio or Shared.Log.Ratio Selection.
#' @param type Two choice for data type(protein and phosphosite). 
#'        Set type = "protein" to handel the protein abundance data;
#'        Set type = "phosphosite" to handel the phosphorylation site abundance data.
#' @return Uniform protein abundance or phosphorylation site abundance data.
#' @author Jun Ge
#' @import stringr
#' @export 
PreprocessCPTAC <- function(exprSet,selection="Unshared",type="protein")
{     
  # Protein abundance data
  if (type == "protein"){
    # Select the shared or unshared data
    if(selection =="Unshared"){
      index <- grepl("Unshared.Log.Ratio",colnames(exprSet))
    } else {
      index <- !grepl("Unshared.Log.Ratio",colnames(exprSet))
    } 
    exprSet <- exprSet[,c(1,which(index))]
    rownames(exprSet) <- exprSet[,"Gene"]
    exprSet[,"Gene"] <- NULL
    exprSet <- exprSet[-c(1:3),]
    # Delete the ".Unshared.Log.Ratio" 
    colnames(exprSet) <- stringr::str_replace(colnames(exprSet)," Log Ratio","") 
    
    # Phosphorylation site abundance data
  } else if (type == "phosphosite"){
    colnames(exprSet) <- stringr::str_replace(colnames(exprSet)," Log Ratio","") 
    exprSet$Phosphosite <- stringr::str_replace(exprSet$Phosphosite,".*:",paste(exprSet$Gene,":",sep = ""))
    exprSet <- exprSet[!duplicated(exprSet$Phosphosite),]
    index <- which(colnames(exprSet) %in% c("Peptide","Gene","Organism"))
    exprSet <- exprSet[,-index]
    rownames(exprSet) <- exprSet[,"Phosphosite"]
    exprSet[,"Phosphosite"] <- NULL
  }
  return(exprSet)
}











