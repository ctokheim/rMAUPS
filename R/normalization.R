#' Normalize Fisher's data
#'
#' @param df.prot A matrix-like object.
#' @param norm method for normalization.
#' @param log2 Boolean, whether perform log2 transformation.
#'
#' @return A matrix.
#' @author Wubing Zhang
#' @export
#'
normalizeProteomeDiscoverer <- function(df.prot, norm = "median", log2 = FALSE){
  # Identify and count reporter ion channels
  channels <- grep("^Abundance.F..(126|127|128|129|130|131)", colnames(df.prot), value=TRUE)
  if (length(channels)==0) {
    channels <- grep("^Abundance.F...(126|127|128|129|130|131)", colnames(df.prot), value=TRUE)
  }

  # Identify contrasts from column names
  temp <- unique(gsub(".*Sample.", "", channels))
  contrasts <- temp[!grepl("DMSO", temp)]

  # Subset data frame
  df <- df.prot[, channels]
  rownames(df) <- df.prot$Accession
  colnames(df) <- gsub("^Abundance.F...|Sample.", "", colnames(df))

  # Remove protein contaminants
  idx1 = toupper(df.prot$Contaminant)!="TRUE"
  # Remove data with unique peptides below threshold
  idx2 = df.prot$Number.of.Unique.Peptides >= 2
  df <- df[idx1&idx2, ]

  # Remove data with NAs or low sum of reporter ion intensities
  idx1 <- rowSums(df, na.rm = TRUE)>0
  df <- df[idx1, ]

  normdf = normalizeMS(df, norm, log2)

  # Normalize and scale data to create relative abundance matrix
  # Remove negative values from previous variables
  return(normdf)
}

#' Normalize Fisher's data
#'
#' @param df A matrix-like object.
#' @param norm method for normalization, such as none, median, medianratio, scale, robustz, quantile, loess.
#' @param log2 Boolean, whether perform log2 transformation.
#'
#' @return A matrix.
#' @author Wubing Zhang
#' @export
#'
normalizeMS <- function(df, norm = "median", log2 = FALSE){
  normdf = df
  if(norm=="medianratio"){
    # Median ratio normalization
    geomeans <- exp(rowMeans(log(normdf)))
    Sizes <- apply(normdf, 2, function(cnts)
      median((cnts/geomeans)[geomeans > 0]))
    normdf = t(t(normdf)/Sizes)
    if(log2) normdf = log2(normdf)
  }else if(norm=="median"){
    if(log2){
      mid = apply(log2(normdf), 2, median, na.rm = TRUE)
      normdf = t(t(log2(normdf)) - mid)
    }else{
      mid = apply(normdf, 2, median, na.rm = TRUE)
      normdf = t(t(normdf) - mid)
    }
  }else if(norm=="scale"){
    if(log2){
      normdf = scale(log2(normdf))
    }else{
      normdf = scale(normdf)
    }
  }else if(norm=="robustz"){
    if(log2){
      M = apply(log2(normdf), 2, median, na.rm = TRUE)
      MAD = apply(log2(normdf), 2, mad, na.rm = TRUE)
      normdf = t((t(log2(normdf)) - M) / MAD)
    }else{
      M = apply(normdf, 2, median, na.rm = TRUE)
      MAD = apply(normdf, 2, mad, na.rm = TRUE)
      normdf = t((t(normdf) - M) / MAD)
    }
  }else if(norm=="quantile"){
    normdf = limma::normalizeQuantiles(normdf)
    if(log2) normdf = log2(normdf)
  }else if(norm=="loess"){
    normdf = limma::normalizeCyclicLoess(normdf)
  }
  return(normdf)
}
