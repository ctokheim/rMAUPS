#' Quality control of proteomic data
#' @param expr The expression profile.
#' @return a ggplot object.
#' @import ggplot2 ggpubr
#' @author Wubing Zhang
#'
#'
countNA <- function(expr
                    # , filename = NULL, width = 4, height = 3, ...
                    ){
  require(ggplot2)
  require(ggpubr)
  if(!class(expr) %in% c("data.frame", "matrix")){
    exprSet <- read.table(expr, sep = "\t", header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  }else{
    exprSet = expr
  }
  ##  Count NA
  genecount <- rowSums(!is.na(exprSet))
  gg = sapply(1:ncol(exprSet), function(x) sum(genecount>x))
  gg = data.frame(Ratio = (1:ncol(exprSet)),
                  Count = gg, stringsAsFactors = FALSE)
  p = ggplot(gg, aes(Ratio, Count))
  p = p + geom_point(color = "gray50")
  p = p + theme_pubr()
  p = p + labs(x = "Number of sample", y = "Number of gene")
  # if(!is.null(filename))
  #   ggsave(filename, p, width = width, height = height, ...)
  return(p)
}
