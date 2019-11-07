#' Filter out rows with NA or low value
#'
#' @param m A matrix-like object.
#' @param minS A numeric, specifying the minimum proportion/number of values should
#' be quantified for each row.
#' @param out A vector, specifying low-quality values (NA and 0).
#' @param impute Imputation method.
#' @return A matrix with the same columns as input matrix.
#' @author Wubing Zhang
#' @export
#'
filterN <- function(m, minS = 3, out = NA, impute = "none"){
  m = as.matrix(m)

  # if(is.null(design)){
  if(minS<1) minS = minS*ncol(m)
  idx = is.na(m)
  if(length(out[!is.na(out)])>0)
    idx = idx | (m==max(out[!is.na(out)]))
  sel = rowSums(!idx)>=minS
  m = imputeNA(m[sel,], method = impute)
  return(m)
}

