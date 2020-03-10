#' Normalization of ProteomeDiscoverer-exported proteomic data
#'
#' @param prot File path (or data frame) of proteomic data,
#' which includes columns of "Contaminant", "Number of unique peptides",
#' "Gene Symbol", and protein abundance of samples.
#' @param norm method for normalization.
#' @param log2 Boolean, whether perform log2 transformation.
#'
#' @return A matrix or output a csv file to local folder.
#' @author Wubing Zhang
#' @export
#'
normalizeProteomeDiscoverer <- function(prot, norm = "median", log2 = FALSE,
                                        outfile = NULL){
  if(length(prot)==1 && dir.exists(prot)){# Process all export data in a folder
    files = list.files(prot, pattern = "export.*txt$", full.names = TRUE)
    for(f in files){
      outfile = basename(gsub("export.*txt$", "normdata.csv", f))
      normalizeProteomeDiscoverer(f, norm, log2, outfile)
    }
  }else if(length(prot)==1 && file.exists(prot)){# Process an export data
    message(format(Sys.time(), "%b-%d-%Y %X "), "Processing ", prot)
    prot = read.table(prot, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, comment.char = "")
  }
  # Remove protein contaminants
  idx1 = toupper(prot$Contaminant)!="TRUE"
  # Remove data with unique peptides below threshold
  idx2 = prot[,grep("Unique.*Peptides",
                       colnames(prot), value = TRUE)]>=2
  prot <- prot[idx1&idx2, ]

  # Identify and count reporter ion channels
  channels <- grep("^Abundance.*Sample", colnames(prot), value=TRUE)
  df <- prot[, channels]

  genename = grep("symbol",colnames(prot), value = TRUE, ignore.case = TRUE)
  dupgenes = prot[, genename][duplicated(prot[, genename])]
  dupdf = t(sapply(unique(dupgenes), function(g){
    idx = prot[, genename]==g
    colMeans(df[idx, , drop = FALSE], na.rm = TRUE)
  }))
  idx = prot[, genename] %in% dupgenes
  df = rbind.data.frame(df[!idx, ], dupdf)
  rownames(df) <- c(prot[!idx, genename], rownames(dupdf))
  colnames(df) <- gsub("Abundance.|Sample.", "", colnames(df))

  # Remove data with NAs or low sum of reporter ion intensities
  idx1 <- rowSums(df, na.rm = TRUE)>0
  df <- df[idx1, ]

  normdf = normalizeMS(df, norm, log2)
  normdf = as.data.frame(normdf, stringsAsFactors = FALSE)
  if(!is.null(outfile)){
    write.csv(normdf, outfile, row.names = TRUE, quote = FALSE)
    return(TRUE)
  }else{
    return(normdf)
  }
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
normalizeProteomics <- function(df, norm = "median", log2 = FALSE){
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
