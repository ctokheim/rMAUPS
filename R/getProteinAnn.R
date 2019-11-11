#' Gene ID conversion between protein ids and gene ids
#'
#' @docType methods
#' @name TransProteinID
#' @rdname TransProteinID
#'
#' @param ids A character vector, input protein identifiers.
#' @param fromType The input ID type, one of "uniprot"(default), "refseq", "enst", "symbol".
#' @param toType The output ID type, similar to `fromType`.
#' @param organism "hsa"(default) or "mmu".
#' @param ensemblHost String, specifying ensembl host, you can use `listEnsemblArchives()`
#' to show all available Ensembl archives hosts.
#' @param update Boolean, specifying whether update built-in protein annotation (needs network and takes time).
#'
#' @return A character vector, named by unique input gene ids.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link[MAGeCKFlute]{TransGeneID}}
#'
#' @examples
#'
#' @import biomaRt
#' @export

TransProteinID <- function(ids, fromType="uniprot", toType="symbol",
                        organism = "hsa", ensemblHost = "www.ensembl.org", update = FALSE){
  requireNamespace("biomaRt")

  #### Verify  parameters ####
  ids = as.character(ids)
  fromType = tolower(fromType)
  toType = tolower(toType)
  if(length(ids)<1) return(c())
  keggcode = rep(c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"), 2)
  names(keggcode) = c(tolower(c("Human", "Mouse", "Rat", "Bovine", "Canine", "Chimp", "Pig")),
                      c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc"))
  if(!tolower(organism)%in%names(keggcode)) stop("Doesn't surport this organism ...")
  # if(!tolower(fromOrg)%in%names(keggcode)) stop("fromOrg error ...")
  # if(!tolower(toOrg)%in%names(keggcode)) stop("toOrg error ...")

  # organism = keggcode[tolower(organism)]
  # fromOrg = keggcode[tolower(fromOrg)]
  # toOrg = keggcode[tolower(toOrg)]

  #### Read annotation file ####
  if(all(c(fromType, toType) %in% c("uniprot", "symbol", "enst", "refseq"))){
    ann <- getProteinAnn(organism, update)
    ann0 = ann; ann0$Uniprot = gsub("-.*", "", ann0$Uniprot)
    ann = rbind.data.frame(ann, ann0)
    if(fromType=="uniprot" | toType=="uniprot"){
      ann = ann[ann$DB %in% c(fromType, toType), c("Uniprot", "ID")]
    }else{
      tmp1 = ann[ann$DB==fromType, ]
      tmp2 = ann[ann$DB==toType, ]
      colnames(tmp1)[3] = fromType
      colnames(tmp2)[3] = toType
      tmp2 = tmp2[!duplicated(tmp2$Uniprot), ]
      ann = merge.data.frame(tmp1[,-2], tmp2[,-2], by = "Uniprot", all = TRUE)
      ann = ann[, -1]
    }
  }else{
    datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                        "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
    ds = datasets[grepl(organism, datasets)]
    ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds, host = ensemblHost)
    ## decide the attributes automatically
    attrs = listAttributes(ensembl)$name
    if(sum(attrs==fromType)==0){
      idx1 = grepl(tolower(fromType), attrs)
      idx = idx1
      if(sum(idx1)>2) idx = idx1&grepl("_id", attrs)
      fromType = ifelse(sum(idx)>0, attrs[idx][1], attrs[idx1][1])
      if(fromType=="hgnc_symbol" & organism=="mmu") fromType = "mgi_symbol"
    }
    if(sum(attrs==toType)==0){
      idx1 = grepl(tolower(toType), attrs)
      idx = idx1
      if(sum(idx1)>2) idx = idx1&grepl("_id", attrs)
      toType = ifelse(sum(idx)>0, attrs[idx][1], attrs[idx1])
      if(toType=="hgnc_symbol" & organism=="mmu") toType = "mgi_symbol"
    }
    ## retrieve the data
    ann = getBM(attributes=c(fromType, toType), mart = ensembl,
                filters = fromType, values = ids)
  }
  ## merge the annotation
  colnames(ann) = c(fromType, toType)
  ## Retain unique conversion
  idx = duplicated(ann[, fromType]) | is.na(ann[, fromType])
  convert = ann[!idx, toType]
  names(convert) = ann[!idx, fromType]
  gene_after = as.character(convert[ids])
  gene_after[is.na(gene_after)] = convert[gsub("\\..*|-.*", "", ids[is.na(gene_after)])]

  names(gene_after) = ids
  return(gene_after)
}

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
                         "&format=tab&force=true&columns=id,entry%20name,reviewed,protein%20names,genes,organism,database(Ensembl),comment(SUBCELLULAR%20LOCATION),database(RefSeq)&sort=score")

  if((!file.exists(locfname)) | update){
    ## Download protein annotation information from uniprot
    download.file(uniprot_link, locfname, quiet = TRUE)
  }

  ## Reorder the mapping file
  uniprot_ann = read.csv(gzfile(locfname), sep = "\t", header = TRUE,
                         quote = "", stringsAsFactors = FALSE, comment.char = "")
  suppressWarnings(try(file.remove(locfname), silent = TRUE))
  colnames(uniprot_ann) = c("Entry", "EntryName", "Status", "Name", "Gene",
                            "Organism", "Ensembl", "Subcellular", "RefSeq")
  summary = apply(uniprot_ann, 1, function(x){
    Uniprot = x[1]
    Symbols = unlist(strsplit(gsub(";$", "", x[5]), " "))
    Symbols[!grepl("\\[", Symbols)] = paste0(Symbols[!grepl("\\[", Symbols)], " [", Uniprot, "]")
    ENSTs = unlist(strsplit(gsub(";$", "", x[7]), ";"))
    ENSTs[!grepl("\\[", ENSTs)] = paste0(ENSTs[!grepl("\\[", ENSTs)], " [", Uniprot, "]")
    Subcellular = gsub("^ |SUBCELLULAR LOCATION: ", "", x[8])
    Subcellular[!grepl("\\[", Subcellular)] = paste0(Subcellular[!grepl("\\[", Subcellular)], " [", Uniprot, "]")
    RefSeq = unlist(strsplit(gsub(";$", "", x[9]), ";"))
    RefSeq[!grepl("\\[", RefSeq)] = paste0(RefSeq[!grepl("\\[", RefSeq)], " [", Uniprot, "]")
    as.vector(rbind(c(Symbols, ENSTs, RefSeq, Subcellular),
          rep(c("symbol", "enst", "refseq", "cellularloc"), c(length(Symbols), length(ENSTs),
                                               length(RefSeq), length(Subcellular)))))
  })
  summary = matrix(unlist(summary), ncol = 2, byrow = TRUE)
  summary = as.data.frame(summary, stringsAsFactors = FALSE)
  colnames(summary) = c("ID", "DB")
  summary$Uniprot = gsub(".*\\[|\\]$", "", summary$ID)
  summary$ID = gsub(" \\[.*", "", summary$ID)
  summary = summary[summary$ID!="", ]

  canonical = readRDS(file.path(system.file("extdata", package = "rMAUPS"),
                                 paste0("uniprot_canonical_", org, ".rds")))
  # canonical = unlist(canonical)
  # names(canonical) = gsub("-.*$", "", canonical)
  # canonical[!grepl("-", canonical)] = paste0(canonical[!grepl("-", canonical)], "-1")
  # saveRDS(canonical, file.path(system.file("extdata", package = "rMAUPS"),
  #                              paste0("uniprot_canonical_", org, ".rds")))
  summary$Uniprot[summary$Uniprot%in%names(canonical)] = canonical[summary$Uniprot[summary$Uniprot%in%names(canonical)]]
  summary = summary[, c("Uniprot", "DB", "ID")]
  saveRDS(summary, rdsann)
  return(summary)
}

