library(biomaRt)
library(stringr)
library(purrr)
library(dplyr)

# parse mupit response
library(httr)
library(jsonlite)

saveSeq <- function() {
  # fetch protein sequence
  ensembl.human = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ann <- getBM(attributes=c('uniprotswissprot'), mart=ensembl.human)
  human.prot = getSequence(id=ann$uniprotswissprot,
                           mart=ensembl.human,
                           seqType=c("peptide"),
                           type="uniprotswissprot")

  # reformat result and calculate protein lengths
  result <- as.data.frame(cbind(human.prot[2], human.prot[1]))
  result['protein_sequence'] <- str_sub(result$peptide, start=1, end=-2)
  result$length <- nchar(result$peptide) - 1

  # write out
  write.table(x=result[,c('uniprotswissprot', 'protein_sequence')], sep='\t', row.names=F, quote=F, file="uniprot_protein_sequence.txt")

  return(result)
}

#myresult <- saveSeq()

# search protein sequences
searchProtSeq <- function(prot_seq_df, regex){
  # regex search all protein sequences
  motif_hits <- str_locate_all(prot_seq_df$protein_sequence, regex) 
  
  # format motif hits
  motif_hits <- map(motif_hits, as.data.frame)
  names(motif_hits) <- prot_seq_df$ID
  motif_hits <- imap(motif_hits, ~.x %>% mutate(UniprotId = .y))
  motif_hits_df <- reduce(compact(motif_hits), rbind)

  return(motif_hits_df)
}


browseProtStructure <- function(protId, start, end, 
                                baseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/?gm=',
                                checkBaseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/rest/showstructure/check?pos='){
  require(httr)
  require(jsonlite)

  # build uniprot query
  uniProtSeq <- ""
  for (pos in seq(start, end)){
    tmp <- paste(protId, pos, sep=':')
    if (uniProtSeq==''){
      uniProtSeq <- tmp
    } else {
      uniProtSeq <- paste(uniProtSeq, tmp, sep=',')
    }
  }

  # check if prot structure available
  jsonResponse <- GET(paste0(checkBaseUrl, uniProtSeq, '&protquery=y'))
  jsonResponseParsed <- content(jsonResponse, as="parsed")

  # construct url
  fullUrl <- paste0(baseUrl, uniProtSeq, '&protquery=y')
  # browse url if there is available prot structure
  if (jsonResponseParsed$hit){
    print(fullUrl)
    browseURL(fullUrl)
  } else {
    print('No protein structure available!') 
  }
}
