#' GetCptacRNAseqData: Together the RNAseq data(Count,FPKM) from NCI Genomic Data Commons (GDC)
#' @description : UCEC,Kidney and Lung cancer RNAseq data are stored in the NCI Genomic Data Commons (GDC).
#'               The data which we get from the GDC are separately.We can use GetCptacRNAseqData
#'               to merged the data into one set and transform the "ensembl_gene_id_version" to "hgnc_symbol".
#' @param gdcdown.path: The path to the GDC data.
#' @param anno_file: Annotation file which is downloaded from GDC.
#' @param type: Data type.E.g."STARCount","HTSeqCount", "HTSeqFPKM", "HTSeqFPKMUQ".
#' @return Merged data, which can do the downstream analysis.
ProcessGDCRNAseq <- function(path,anno_file,type="STARCount"){
  res <- list()
  dirs <- list.files(path)
  dirs <- dirs[!grepl(".txt",dirs)] # remove the "MANIFEST.txt" file
  # get the expression file
  rownames(anno_file) <- anno_file[,"File ID"]
  anno_file[,"File ID"] <- NULL
  anno_file <- anno_file[dirs,]
  anno_file$Type <- "N"
  anno_file$Type[grepl("Tumor",anno_file$`Sample Type`)] <- "T"
  anno_file$`Case ID` <- paste(anno_file$`Case ID`, anno_file$Type, sep = ".")
  anno_file$`Case ID` <- gsub(".*, ","",anno_file$`Case ID`)
  files <- paste(path,dirs,anno_file[,"File Name"],sep = "/")
  if (type=="STARCount"){
    df <- lapply(files,function(x) read.csv(file=x, sep = "\t",stringsAsFactors = F,header=T)) %>% bind_cols()
    df <- df[-c(1:4),]
    gene_expr_which <- c(1,which(grepl("unstranded",colnames(df))))
    gene_expr <- df[,gene_expr_which]
  }else{ # HTSeqCount; HTSeqFPKM; HTSeqFPKMUQ
    df <- lapply(files,function(x) read.csv(file=x, sep = "\t",stringsAsFactors = F,header=F)) %>% bind_cols()
    gene_expr_which <- c(1,seq(2,ncol(df),2))
    gene_expr <- df[,gene_expr_which]
  }
  colnames(gene_expr) <- c("entrze",rownames(anno_file))
  colnames(gene_expr) <- c("entrze",anno_file[,"Case ID"])
  gene_expr$entrze <- str_replace_all(gene_expr$entrze,"\\..*","")
  gene_expr <- gene_expr[,!duplicated(colnames(gene_expr))]
  # Gene ID transform.
  ense_sym <- MAGeCKFlute::TransGeneID(gene_expr$entrze, fromType = "Ensembl", toType = "Symbol",
                          organism = "hsa", useBiomart = FALSE,
                          ensemblHost = "www.ensembl.org")
  ense_sym <- as.data.frame(ense_sym)
  ense_sym <- na.omit(ense_sym)
  #
  gene_expr_symbol <- merge(gene_expr,ense_sym,by.x=1,by.y=0)
  rownames(gene_expr_symbol) <- gene_expr_symbol[,'ense_sym']
  gene_expr_symbol[,'ense_sym'] <- NULL
  gene_expr_symbol[,'entrze'] <- NULL
  gene_expr_symbol <- gene_expr_symbol[order(rownames(gene_expr_symbol)),]
  #================
  normal <- gene_expr_symbol[,grep("\\.N",colnames(gene_expr_symbol))]
  colnames(normal) <- gsub("\\.N","",colnames(normal))
  tumor <- gene_expr_symbol[,grep("\\.T",colnames(gene_expr_symbol))]
  colnames(tumor) <- gsub("\\.T","",colnames(tumor))
  res[[1]] <- normal
  res[[2]] <- tumor
  return(res)
}
