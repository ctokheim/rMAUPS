library(biomaRt)
library(stringr)
library(purrr)
library(dplyr)

# parse mupit response
library(httr)
library(jsonlite)

# search protein sequences
searchProtSeq <- function(prot_seq_df, regex){
  # regex search all protein sequences
  motif_hits <- str_locate_all(prot_seq_df$protein_sequence, regex) 
  
  # format motif hits
  motif_hits <- map(motif_hits, as.data.frame)
  names(motif_hits) <- prot_seq_df$ID
  motif_hits <- imap(motif_hits, ~.x %>% mutate(UniprotId = .y))
  motif_hits_df <- reduce(compact(motif_hits), rbind)

  myorder <- c('UniprotId', 'start', 'end')
  return(motif_hits_df[,myorder])
}


browseProtStructure <- function(protId, start, end, 
                                doBrowse=TRUE,
                                baseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/?gm=',
                                checkBaseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/rest/showstructure/check?pos='){
  requireNamespace(httr)
  requireNamespace(jsonlite)

  stopifnot(length(start)==length(end))

  uniProtSeq <- ""
  for (i in 1:length(start)){
    # build uniprot query
    s <- start[i] ; e <- end[i]
    for (pos in seq(s, e)){
      tmp <- paste(protId, pos, sep=':')
      if (uniProtSeq==''){
        uniProtSeq <- tmp
      } else {
        uniProtSeq <- paste(uniProtSeq, tmp, sep=',')
      }
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
    if (!doBrowse) {
      # pass, do nothing 
    } else if (Sys.getenv('R_BROWSER')!="") {
      browseURL(fullUrl)
    } else {
      print('Browser not set!')
    }
  } else {
    print('No protein structure available!') 
  }
}
