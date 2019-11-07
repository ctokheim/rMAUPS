
#' Retrieve gene annotations from the NCBI, HNSC, and Uniprot databases.
#'
#' @docType methods
#' @name getProteinAnn
#' @rdname getProteinAnn
#'
#' @param org Character, hsa (default), bta, cfa, mmu, ptr, rno, ssc are optional.
#' @param update Boolean, indicating whether download current annotation.
#' @return A data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' ann = getProteinAnn("hsa")
#' head(ann)
#'
#' @export
#'
getProteinAnn <- function(org = "hsa", update = FALSE){
  #### Read rds file directly ####
  rdsann = file.path(system.file("extdata", package = "rMAUPS"),
                     paste0("Proteome_Annotation_", org, ".rds"))
  if(file.exists(rdsann) & !update) return(readRDS(rdsann))

  #### NCBI gene annotation ####
  proteome_code = c("hsa"="up000005640", "mmu"="up000000589")
  locfname <- file.path(system.file("extdata", package = "rMAUPS"),
                        paste0("uniprot_proteome_", proteome_code[org], ".tab"))
  uniprot_link <- paste0("https://www.uniprot.org/uniprot/?query=proteome:", proteome_code[org],
                         "&format=tab&force=true&columns=id,entry%20name,reviewed,protein%20names,genes,organism,database(Ensembl),comment(SUBCELLULAR%20LOCATION)&sort=score")

  if((!file.exists(locfname)) | update){
    ## Download protein annotation information from uniprot
    download.file(uniprot_link, locfname, quiet = TRUE)
  }

  ## Reorder the mapping file
  ncbi_ann = read.csv(gzfile(locfname), sep = "\t", header = TRUE,
                      quote = "", stringsAsFactors = FALSE, comment.char = "")
  ncbi_ann = ncbi_ann[, c("GeneID", "Symbol", "Synonyms", "dbXrefs", "type_of_gene", "description")]
  colnames(ncbi_ann)[c(1,2,6)] = c("entrez", "symbol", "fullname")
  ncbi_ann$hgnc = gsub("\\|.*", "", gsub(".*HGNC:", "", ncbi_ann$dbXrefs))
  ncbi_ann$ensembl = gsub("\\|.*", "", gsub(".*Ensembl:", "", ncbi_ann$dbXrefs))
  ncbi_ann$hgnc[!grepl("HGNC", ncbi_ann$dbXrefs)] = ""
  ncbi_ann$ensembl[!grepl("Ensembl", ncbi_ann$dbXrefs)] = ""

  synonyms_row = matrix(unlist(apply(ncbi_ann, 1, function(x){
    tmp = unlist(strsplit(x[3], "[|]"))
    if(tmp[1]!="" & tmp[1]!="-") return(as.vector(rbind(x[1], tmp, x[7], x[8], x[6])))
    return(NULL)
  })) , ncol=5, byrow = TRUE)
  colnames(synonyms_row) = c("entrez", "symbol", "hgnc", "ensembl", "fullname")
  ncbi_ann = rbind(ncbi_ann[,c(1,2,7,8,6)], synonyms_row)
  ncbi_ann = ncbi_ann[,-5]

  #### HGNC gene annotation ####
  if(org=="hsa"){
    locfname2 = file.path(system.file("extdata", package = "MAGeCKFlute"), "HGNC_GeneID_annotation.txt.gz")
    if((!file.exists(locfname2)) | update){
      ## Download gene information from HGNC
      refname <- "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_aliases&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&col=gd_pub_eg_id&col=md_refseq_id&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
      download.file(refname, locfname2, quiet = TRUE)
    }
    ## Reorder the mapping file
    hgnc_ann = read.csv(gzfile(locfname2), sep = "\t", header = TRUE,
                        quote = "", stringsAsFactors = FALSE, comment.char = "")
    hgnc_ann = hgnc_ann[, c(8,2,4,1,7,3)]
    names(hgnc_ann) = c("entrez", "symbol", "synonyms", "hgnc", "ensembl", "fullname")
    hgnc_ann$hgnc = gsub("HGNC:", "", hgnc_ann$hgnc)
    synonyms_row = matrix(unlist(apply(hgnc_ann, 1, function(x){
      tmp = unlist(strsplit(x[3], ", "))
      if(length(tmp)>0) return(as.vector(rbind(x[1], tmp, x[4], x[5], x[6])))
      return(NULL)
    })) , ncol=5, byrow = TRUE)
    colnames(synonyms_row) = c("entrez", "symbol", "hgnc", "ensembl", "fullname")
    hgnc_ann = rbind(hgnc_ann[,-3], synonyms_row)
    hgnc_ann = hgnc_ann[, -5]
  }

  #### Ensembl gene annotation ####
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  ds = datasets[grepl(org, datasets)]
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
  symbol <- ifelse(org=="mmu", "mgi_symbol", "hgnc_symbol")
  ensembl_ann = getBM(attributes=c("entrezgene_id", symbol, "hgnc_id", "ensembl_gene_id"), mart = ensembl)
  colnames(ensembl_ann) = c("entrez", "symbol", "hgnc", "ensembl")
  ensembl_ann$hgnc = gsub("HGNC:", "", ensembl_ann$hgnc)

  #### Merge HGNC and NCBI annotation ####
  if(org=="hsa"){
    data = rbind.data.frame(ensembl_ann, hgnc_ann, ncbi_ann)
  }else data = rbind.data.frame(ensembl_ann, ncbi_ann)
  idx = duplicated(paste(data$entrez, data$symbol, data$ensembl, sep = "_"))
  data = data[!idx, ]
  data$entrez = as.integer(data$entrez)
  data$hgnc = as.integer(data$hgnc)
  data = data[order(data$entrez), ]
  rownames(data) = NULL
  saveRDS(data, rdsann)
  return(data)
}

