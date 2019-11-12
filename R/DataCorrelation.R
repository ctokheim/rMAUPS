#' Correlation between RNAseq and Protein
#' @docType methods
#' @name DataCorrelation
#' @rdname DataCorrelation
#' @param ds1 Path to protein abundance data
#' @param ds2 Path to RNAseq data
#' @param method One of "pearson", "kendall", or "spearman" (default).
#' @return A ggplot instance.
#' @import ggplot2 ggpubr
#' @author Wubing Zhang

DataCorrelation <- function(ds1, ds2, method="spearman"
                            # , filename = NULL, width = 5, height = 4, ...
                            ){
  if(!class(ds1) %in% c("data.frame", "matrix")){
    Protein <- read.table(ds1, sep = "\t",stringsAsFactors = FALSE,
                          check.names = FALSE, quote = "", row.names = 1)
  }else{
    Protein = ds1
  }
  Protein <- na.omit(Protein)
  if(!class(ds2) %in% c("data.frame", "matrix")){
    RNA <- read.table(ds2, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, quote = "", row.names = 1)
  }else{
    RNA = ds2
  }
  RNA <- na.omit(RNA)

  overlap_gene <- intersect(rownames(Protein), rownames(RNA))
  com_gene <- length(overlap_gene)
  overlap_sample <- intersect(colnames(Protein), colnames(RNA))
  com_sample <- length(overlap_sample)
  if(com_gene<100 | com_sample<4){
    message("Number of overlap genes: ", com_gene,
            "\nNumber of overlap samples: ", com_sample)
    stop("limited number of samples or genes!")
  }
  Protein_sub <- Protein[overlap_gene, overlap_sample]
  RNA_sub <- RNA[overlap_gene, overlap_sample]

  ## Sample correlation
  corvalues_s = diag(cor(Protein_test, RNA_test, method = method))
  corvalues_g = diag(cor(t(Protein_test), t(RNA_test), method = method))
  gg = data.frame(Cor = c(corvalues_s, corvalues_g),
                  Type = rep(c("Sample", "Gene"), c(com_sample, com_gene)))
  p = ggplot(gg, aes(Cor, color = Type))
  p = p + geom_density()
  p = p + labs(x = "Correlation", y = "Density")
  p = p + theme_pubr()
  # if(!is.null(filename))
  #   ggsave(filename, p, width = width, height = height, ...)
  return(p)
}
