#' Differential expression analysis
#'
#' @docType methods
#' @name DEAnalyze
#' @rdname DEAnalyze
#'
#' @param obj Matrix like object or an ExpressionSet instance.
#' @param SampleAnn Matrix like object (only when obj is a matrix),
#' the rownames should match colnames in obj, and the first column should be Condition.
#' @param type "Array", "RNASeq" or "msms", only needed when obj is matrix like object.
#' @param method Differential expression analysis method, e.g. limma, DESeq2, GFOLD,
#' glm.pois, glm.qlll, and glm.nb.
#' @param paired Boolean, specifying whether perform paired comparison.
#' @param GeneAnn Matrix like object (only when obj is a matrix), the rownames should match rownames of obj.
#' @param return Character, either data.frame or ExpressionSet, specifying the return object type.
#' @param app.dir The path to application (e.g. GFOLD).
#'
#' @return An ExpressionSet instance.
#' @seealso \code{\link{ExpressionSet-class}}
#'
#' @author Wubing Zhang
#'
#' @import Biobase
#' @importFrom limma lmFit eBayes topTable voom
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq lfcShrink
#' @importFrom msmsTests msms.glm.pois msms.glm.qlll msms.edgeR
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom MSnbase MSnSet
#' @importFrom msmsEDA pp.msms.data
#' @import utils stats
#' @export

DEAnalyze <- function(obj, SampleAnn = NULL, type = "Array",
                      method = "limma", paired = FALSE, GeneAnn = NULL,
                      return = c("data.frame", "ExpressionSet")[1],
                      app.dir = "/Users/Wubing/Applications/gfold/gfold"){
  #### Create a new object ####
  if(is.matrix(obj) | is.data.frame(obj)){
    obj = stats::na.omit(obj)
    colnames(SampleAnn)[1] = "Condition"
    if(paired) colnames(SampleAnn)[2] = "Sibs"
    samples = intersect(rownames(SampleAnn), colnames(obj))
    if(length(samples)<2) stop("Too small sample size !!!")
    expr <- as.matrix(obj[, samples])
    SampleAnn = SampleAnn[samples, , drop = FALSE]
    # obj = new("ExprDataSet", rawdata = expr, SampleAnn = SampleAnn, type = type)
    obj = ExpressionSet(assayData = expr, phenoData = AnnotatedDataFrame(SampleAnn))
    if(!is.null(GeneAnn)){
      tmp = GeneAnn[rownames(obj), drop = FALSE]
      rownames(tmp) = rownames(obj)
      slot(obj, "featureData") = AnnotatedDataFrame(tmp)
    }
  }

  #### Build design matrix ####
  if(paired){
    Sibs = factor(pData(obj)$Sibs)
    Condition = factor(pData(obj)$Condition)
    design = stats::model.matrix(~Sibs+Condition)
  }else{
    design = stats::model.matrix(~1+Condition, pData(obj))
    rownames(design) = sampleNames(obj)
  }
  if(tolower(type) == "array"){
    requireNamespace("limma")
    # exprs(obj) = normalizeQuantiles(exprs(obj))
    #"ls" for least squares or "robust" for robust regression
    fit = eBayes(lmFit(exprs(obj), design, na.rm=TRUE))
    res = topTable(fit, adjust.method="BH", coef=ncol(design), number = Inf)
    res = res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
    colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
  }else if(tolower(type) == "rnaseq"){
    if(tolower(method) == "deseq2"){
      requireNamespace("DESeq2")
      # exprs(obj) = TransformCount(exprs(obj), method = "vst")
      # DESeq2
      dds = DESeqDataSetFromMatrix(exprs(obj), colData = pData(obj), design = design)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::lfcShrink(dds, coef = ncol(design), quiet = TRUE)
      res$padj[is.na(res$padj)] = 1
      res = res[, c("log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(tolower(method) == "limma"){
      requireNamespace("limma")
      # exprs(obj) = TransformCount(exprs(obj), method = "voom")
      # limma:voom
      dge <- DGEList(counts=exprs(obj))
      dge <- calcNormFactors(dge)
      dge <- voom(dge, design, plot=FALSE)
      fit <- eBayes(lmFit(dge, design))
      res = topTable(fit, adjust.method="BH", coef=ncol(design), number = nrow(exprs(obj)))
      res = res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(tolower(method) == "edger"){
      # exprs(obj) = TransformCount(exprs(obj), method = "voom")
      dge <- edgeR::DGEList(counts=exprs(obj))
      dge <- edgeR::calcNormFactors(dge)
      dge <- edgeR::estimateDisp(dge, design, robust=TRUE)
      fit <- edgeR::glmFit(dge, design)
      lrt <- edgeR::glmLRT(fit)
      res <- edgeR::topTags(lrt, n = nrow(exprs(obj)))
      res = res$table[, c("logFC", "logCPM", "logFC", "PValue", "FDR")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(tolower(method) == "gfold"){
      exprs(obj) = TransformCount(exprs(obj), method = "voom")
      # GFOLD
      tmp = mapply(function(x){
        write.table(cbind(NA, exprs(obj)[,x], NA, NA),
                    file=paste0(colnames(exprs(obj))[x], ".txt"),
                    sep="\t", col.names=FALSE)}, x=1:ncol(exprs(obj)))
      lev = levels(pData(obj)$Condition)
      ctrlname = rownames(pData(obj))[pData(obj)$Condition==lev[1]]
      treatname = rownames(pData(obj))[pData(obj)$Condition==lev[2]]
      system(paste0(app.dir, " diff -s1 ", paste0(ctrlname, collapse=","),
                    " -s2 ", paste0(treatname, collapse=","), " -suf .txt -o gfold_tmp")
             )
      res = read.table("gfold_tmp", row.names=1, stringsAsFactors = FALSE)
      res = res[, c(4, 5, 2, 3, 3)]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
      tmp = file.remove(paste0(ctrlname, ".txt"), paste0(treatname, ".txt"),
                        "gfold_tmp", "gfold_tmp.ext")
    }else{
      stop("Method not available for RNA-seq data !!!")
    }
  }else if(tolower(type) == "msms"){
    if (tolower(method) == "limma"){
      requireNamespace("limma")
      #"ls" for least squares or "robust" for robust regression
      fit = eBayes(lmFit(exprs(obj), design))
      res = topTable(fit, adjust.method="BH", coef=ncol(design), number = Inf)
      res = res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(grepl("^glm\\.", tolower(method))){
      requireNamespace("msmsTests")
      fd <- data.frame(gene = rownames(exprs(obj)),
                       row.names = rownames(exprs(obj)),
                       stringsAsFactors = FALSE)
      MSnSet_obj <- MSnSet(exprs=exprs(obj), fData=fd,
                           pData=pData(obj))
      MSnSet_obj <- pp.msms.data(MSnSet_obj)  # pp.msms.data function used to deleted genes which all expression is 0.

      null.f <- "y~1"
      alt.f <- "y~Condition"
      div <- colSums(exprs(obj), na.rm = TRUE)
      ### msmsTests method
      if(tolower(method)=="glm.pois"){
        res <- msms.glm.pois(MSnSet_obj, alt.f, null.f, div=div)
      }else if(tolower(method)=="glm.qlll"){
        res <- msms.glm.qlll(MSnSet_obj, alt.f, null.f, div=div)
      }else if(tolower(method)=="glm.nb"){
        res <- msms.edgeR(MSnSet_obj, alt.f, null.f, div=div)
      }
      res$baseMean = rowMeans(exprs(obj))[rownames(res)]
      res$padj = p.adjust(res$p.value, method = "BH")
      res = res[, c(1,4,2:3,5)]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }
  }else{
    stop("Data type error! ")
  }
  res = as.data.frame(res)
  if(return == "data.frame") return(res)
  tmp = cbind(obj@featureData@data, res[rownames(obj@featureData@data), ])
  rownames(tmp) = rownames(obj@featureData@data)
  slot(obj, "featureData") = AnnotatedDataFrame(tmp)
  return(obj)
}

