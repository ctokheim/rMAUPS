#' Principal component visualization
#' @param mat A data matrix.
#' @param color The column name specifying the color.
#' @param label The column name specifying the label.
#' @param filename The file name of output figure.
#' @param width The width of the output figure.
#' @param height The height of the output figure.
#' @param ... parameters in ggsave.
#'
#' @return a ggplot instance.
#' @import ggplot2 ggpubr
#' @author Wubing Zhang
#' @export
#'
pcView <- function(mat, color = gsub(".*_", "", colnames(mat)),
                   filename = NULL, width = 5, height = 4, ...){
  require(ggplot2)
  require(ggpubr)
  mat = as.matrix(mat)
  tmp = prcomp(t(mat))
  gg = as.data.frame(tmp$x[,1:2], stringsAsFactors = FALSE)
  gg$Color = color
  p = ggplot(gg, aes(PC1, PC2, color = Color))
  p = p + geom_point()
  p = p + theme_pubr()
  p = p + labs(color = NULL)
  if(!is.null(filename))
    ggsave(filename, p, width = width, height = height, ...)
  return(p)
}
