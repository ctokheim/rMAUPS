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
#'
pcView <- function(mat, color = gsub(".*_", "", colnames(mat)), label = NULL,
                   filename = NULL, width = 5, height = 4, ...){
  require(ggplot2)
  require(ggpubr)

  tmp = prcomp(t(data))
  gg = as.data.frame(tmp$x[,1:2], stringsAsFactors = FALSE)
  if(!(is.null(label)|label=="")){
    if(label==0) gg$Label = rownames(mat)
    else gg$Label = mat[, label]
  }
  if(!(is.null(color)|color=="")){
    if(color==0) gg$Color = rownames(mat)
    else gg$Color = mat[, color]
  }
  p = ggplot(gg, aes(PC1, PC2))
  if(all(c("Label", "Color") %in% colnames(gg))){
    p = ggplot(gg, aes(PC1, PC2, label = Label, color = Color))
  }else if("Color" %in% colnames(gg)){
    p = ggplot(gg, aes(PC1, PC2, color = Color))
  }else if("Label" %in% colnames(gg)){
    p = ggplot(gg, aes(PC1, PC2, label = Label))
  }
  p = p + geom_point()
  if("Label" %in% colnames(gg))
    p = p + ggrepel::geom_text_repel()
  p = p + theme_pubr()
  p = p + labs(color = NULL)
  ggsave(filename, p, width = width, height = height, ...)
  return(p)
}
