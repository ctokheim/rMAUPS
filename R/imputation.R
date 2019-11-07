#' Filter out rows with NA or low value
#'
#' @param m A matrix-like object.
#' @param method method for imputation, such as knn, lowAbundanceResampling, ReplicateBasedResampling
#' @param k Integer, parameter for knn.
#' @param rowmax parameter for knn.
#' @param colmax parameter for knn.
#'
#' @return A matrix.
#' @author Wubing Zhang
#' @export
#'
imputeNA <- function(m, method = "knn", k = 30, rowmax = 0.95, colmax = 0.95, ...){
  require(impute)
  imputed_m = m
  if(tolower(method) == "knn"){
    imputed_m = impute.knn(m, k = k, rowmax = rowmax, colmax = colmax,
                           maxp = floor(nrow(m)/1000)*1000)$data
  }else if(tolower(method) == "lowabundanceresampling"){
    imputed_m = lowAbundanceResampling(m)
  }else if(tolower(method) == "replicatebasedresampling"){
    imputed_m = ReplicateBasedResampling(m)
  }
  return(imputed_m)
}

#' low Abundance Resampling method from protein discover
#' @param df Matrix-like object.
#' @param percent cutoff for low abundance values.
lowAbundanceResampling <- function(df, percent = 0.05){
  distr = df[!is.na(df)]
  df[is.na(df)] = rnorm(sum(is.na(df)), mean(distr), sd(distr))
  return(df)
}

#' Replicate based resampling method from protein discover
#' @param df Matrix-like object.
ReplicateBasedResampling <- function(df){
  replicate = rep(1, ncol(df))
  for(i in unique(replicate)){
    mid = rowMeans(df[, replicate==i], na.rm = TRUE)
    sd = matrixStats::rowSds(df[, replicate==i])
    mod = lm(sd~mid, data.frame(sd = sd[!is.na(sd)], mid = mid[!is.na(sd)]))
    pred = predict(mod, data.frame(mid = mid))
    for(j in which(is.na(sd)&(!is.na(mid)))){
      tmpj = is.na(df[j, replicate==i])
      df[j, which(replicate==i)[tmpj]] = rnorm(sum(tmpj), mid[j], pred[j])
    }
    return(lowAbundanceResampling(df))
  }
}
