#' RabitView
#' @docType methods
#' @name RabitView
#' @rdname RabitView
#' @param rabitRes Rabit Result.
#' @param top Numbers of factors shown in the plot.
#' @return Bar plot
#' @import ggplot2
#' @author Jun Ge
#' @export
RabitView <- function(rabitRes, top=20){
  colnames(rabitRes) <- c("TF","Cp","Estimate","SE","t.value","P.value")
  rabitRes <- rabitRes[-c(1:3),]
  rabitRes$TF <- gsub(".*\\.","",rabitRes$TF)
  rabitRes$TF <- gsub("@.*","",rabitRes$TF)
  rabitRes <- rabitRes[order(rabitRes$P.value,decreasing = FALSE),]
  rabitRes <- rabitRes[!duplicated(rabitRes$TF),]
  rabitRes$logP <- -log10(rabitRes$P.value)
  rabitRes <- rabitRes[1:top, colnames(rabitRes) %in%  c("TF","logP")]
  rabitRes$TF <- factor(rabitRes$TF,levels = rev(rabitRes$TF))
  p <- ggplot(rabitRes, aes(x=TF,y=logP)) + geom_histogram(stat = "identity",fill="black")+theme_bw()
  p <- p + scale_y_continuous(position = "right")
  p <- p + xlab("")+ylab("-Log10(p-value)")
  p <- p+coord_flip()
  p <- p + theme(axis.text.x = element_text(size =10))
  p <- p + theme(axis.text.y = element_text(size = 10))
  p <- p + theme(axis.title.y = element_text(size = 15,face = "bold"))
  p <- p + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"))
  return(p)
}






