#' QC for Batch effect
#' @docType methods
#' @name BatchQC
#' @rdname BatchQC
#' @param exprSet Protein expression dataset.
#' @param batch Samples batch information.(Three column: Channel, Plex, Batch)
#' the rownames should match colnames in exprSet.
#' @return A list containing three pca plot.
#' @import ggplot2
#' @author Jun Ge
#' @export
BatchQC <- function(expeSet, batch){
    common <- intersect(colnames(expeSet),rownames(batch))
    expeSet <- expeSet[,colnames(expeSet) %in% common]
    batch <- batch[rownames(batch) %in% common,]
    expr_pca <- prcomp(t(expeSet))
    pca_out <- as.data.frame(expr_pca$x)
    pca_out$Channel <- as.factor(batch[rownames(pca_out),"Channel"])
    pca_out$Plex <- as.factor(batch[rownames(pca_out),"Plex"])
    pca_out$Batch <- as.factor(batch[rownames(pca_out),"Batch"])
    xlab <- paste("PC1","(",round((summary(expr_pca))$importance[2,1]*100,1),"%)",sep="")
    ylab <- paste("PC2","(",round((summary(expr_pca))$importance[2,2]*100,1),"%)",sep="")
    doPlot <- function(select){
       p <- ggplot(pca_out,aes_string(x="PC1",y="PC2",color=select))+theme_bw()
       p <- p+geom_point()+ xlab(xlab) + ylab(ylab)
       p <- p + theme(axis.line = element_line(size=0.5, colour = "black"),
                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                     panel.border = element_blank(), panel.background = element_blank())
       p <- p + theme(axis.text.x = element_text(size = 10, face = "bold"))
       p <- p + theme(axis.text.y = element_text(size = 10,face = "bold"))

       return(p)
    }
    res <- lapply(c("Channel","Plex","Batch"), doPlot) #
    return(res)
}
