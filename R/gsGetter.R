#' Extract pathway annotation from GMT file.
#'
#' @docType methods
#' @name gsGetter
#' @rdname gsGetter
#'
#' @param gmtpath The path to customized gmt file.
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (KEGG, REACTOME, C2_CP), GO (GOBP, GOCC, GOMF),
#' MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA), C2_CGP),
#' C3 (C3_MIR, C3_TFT), C4, C6, C7, HALLMARK)
#' and Complex (CORUM). Any combination of them are also accessible
#' (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets to load.
#' @param organism 'hsa' or 'mmu'.
#'
#' @return A three-column data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' gene2path = gsGetter(type = "REACTOME+CORUM")
#' head(gene2path)
#'
#' @export
#'
gsGetter <- function(gmtpath = NULL, type = "All", limit = c(0, Inf),
                     organism = 'hsa'){
  ## Normalize type
  type = toupper(unlist(strsplit(type, "\\+")))
  if("ALL" %in% type) type = c("PATHWAY", "GO", "COMPLEX", "MSIGDB")
  if("MSIGDB" %in% type) type = c("C1", "C2", "C3", "C4", "GO", "C6", "C7", "HALLMARK", type)
  if("C2" %in% type) type = c("KEGG", "REACTOME", "C2", type)
  if("GO" %in% type) type = c("GOBP", "GOCC", "GOMF", type)
  if("PATHWAY" %in% type) type = c("KEGG", "REACTOME", "C2_CP", type)
  if("COMPLEX" %in% type) type = c("CORUM", type)
  type = setdiff(type, c("ALL", "MSIGDB", "C2", "GO", "PATHWAY", "COMPLEX"))
  ## read GMT files
  if(!is.null(gmtpath)){
    gene2path = ReadGMT(gmtpath, limit = limit)
  }else{
    gene2path = data.frame()
    gs_prefix = c("go.all.entrez.", "kegg.all.entrez.", "corum.all.entrez.",
                  "reactome.all.entrez.", "msigdb.all.entrez.")
    if(any(c("GOBP", "GOCC", "GOMF") %in% type)){
      gsfile = file.path(system.file("extdata", package = "rMAUPS"),
                         paste0("go.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      tmp = tmp[grepl(paste(paste0("^",type),collapse="|"), tmp$PathwayID), ]
      gene2path = rbind(gene2path, tmp)
    }
    if("KEGG" %in% type){
      gsfile = file.path(system.file("extdata", package = "rMAUPS"),
                         paste0("kegg.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if("REACTOME" %in% type){
      gsfile = file.path(system.file("extdata", package = "rMAUPS"),
                         paste0("reactome.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if("CORUM" %in% type){
      gsfile = file.path(system.file("extdata", package = "rMAUPS"),
                         paste0("corum.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if(length(setdiff(type, c("GOBP", "GOCC", "GOMF", "KEGG", "CORUM", "REACTOME")))>0){
      gsfile = file.path(system.file("extdata", package = "rMAUPS"),
                         paste0("msigdb.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      tmp = tmp[grepl(paste(paste0("^",type),collapse="|"), tmp$PathwayID), ]
      gene2path = rbind(gene2path, tmp)
    }
    gene2path = na.omit(gene2path)
    count_gene = table(gene2path$PathwayID)
    pathways = names(count_gene)[count_gene>limit[1] & count_gene<=limit[2]]
    gene2path = gene2path[gene2path$PathwayID%in%pathways, ]
  }
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  # gene2path$PathwayName = tolower(gsub("_", " ", gene2path$PathwayName))
  gene2path$Gene = as.character(gene2path$Gene)
  return(gene2path)
}
