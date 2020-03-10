#' Quality control of proteomic data
#'
#' @param data A data frame of the proteomics data.
#' @param meta A vector of conditions corresponding to the samples in the proteomics data.
#' @param proj.name A character as the prefix of output files.
#' @param outdir The path to a local directory.
#'
#' @return a list of ggplot objects.
#' @import ggplot2 MAGeCKFlute
#' @author Wubing Zhang
#' @export
#'
qcProteomics <- function(data, condition = NULL, proj.name = NA, outdir = NULL){
  ## QCs associated with missing values
  missflag = FALSE
  if(sum(is.na(data))>0){
    missflag = TRUE
    p3 = countNA(data)
    gg = data.frame(gene = rownames(data), NAs = rowSums(is.na(data)))
    p4 = DensityView(gg[,2,drop=FALSE], xlab = "The number of missing value")
    p4 = p4 + theme(legend.position = "none")
    gg = data.frame(sample = colnames(data), Detection = colSums(!is.na(data)))
    p5 = DensityView(gg[,2,drop=FALSE], xlab = "The number of detected gene")
    p5 = p5 + theme(legend.position = "none")
    p6 = BarView(gg, "sample", "Detection", fill = "#8da0cb",
                 ylab = "The number of detected gene")
    p6 = p6 + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))

    if(!is.null(outdir)){
      ggsave(paste0(outdir, "/", proj.name, "_count_NA.png"), p3, width = 4, height = 3.5)
      ggsave(paste0(outdir, "/", proj.name, "_density_NA_acrossGene.png"), p4, width = 4, height = 3.5)
      ggsave(paste0(outdir, "/", proj.name, "_density_detection.png"), p5, width = 4, height = 3.5)
      ggsave(paste0(outdir, "/", proj.name, "_bar_detection.png"), p6, width = 4, height = 3.5)
    }
  }else{
    p3 = p4 = p5 = p6 = NULL
  }
  ## Consistency between samples
  p1 = ViolinView(data, ylab = "Logarithmic abundance")
  p1 = p1 + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
  p2 = pcView(data, color = condition)
  if(!is.null(outdir)){
    ggsave(paste0(outdir, "/", proj.name, "_violin_abundance.png"), p1, width = 4, height = 3.5)
    ggsave(paste0(outdir, "/", proj.name, "_pcview_samples.png"), p2, width = 4, height = 3.5)
  }
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
}
