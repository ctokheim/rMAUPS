#library(stringr)
#library(purrr)
#library(dplyr)
# parse mupit response
#library(httr)
#library(jsonlite)

#' Search protein sequences using a regular expression
#'
#' @docType methods
#' @name searchProtSeq
#' @rdname searchProtSeq
#'
#' @param prot_seq_df A data.frame containing uniprot IDs under the "ID" column and protein sequence under the "protein_sequence" column
#' @param regex a string containing the regular expression to search against the protein sequences
#'
#' @return A data frame object
#'
#' @author Collin Tokheim
#'
#' @import purrr
#' @import dplyr
#' @import stringr
#' @export
searchProtSeq <- function(prot_seq_df, regex){
  # regex search all protein sequences
  motif_hits <- stringr::str_locate_all(prot_seq_df$protein_sequence, regex) 
  
  # format motif hits
  motif_hits <- purrr::map(motif_hits, as.data.frame)
  names(motif_hits) <- prot_seq_df$ID
  motif_hits <- purrr::imap(motif_hits, ~.x %>% mutate(UniprotId = .y))
  motif_hits_df <- dplyr::reduce(compact(motif_hits), rbind)

  myorder <- c('UniprotId', 'start', 'end')
  return(motif_hits_df[,myorder])
}


#' View protein sequence positions on protein structure
#'
#' @docType methods
#' @name browseProtStructure
#' @rdname browseProtStructure
#'
#' @param protId a uniprot id string
#' @param start a vector of start positions
#' @param end a vector of end positions
#' @param doBrowse a boolean indicator on whether to open browser to view protein structure
#' @param baseUrl string for location of mupit service
#' @param checkBaseUrl string for location of mupit service to check availability of protein structure
#'
#' @return None
#'
#' @author Collin Tokheim
#'
#' @import httr
#' @export
browseProtStructure <- function(protId, start, end, 
                                doBrowse=TRUE,
                                baseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/?gm=',
                                checkBaseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/rest/showstructure/check?pos='){

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
  jsonResponse <- httr::GET(paste0(checkBaseUrl, uniProtSeq, '&protquery=y'))
  jsonResponseParsed <- httr::content(jsonResponse, as="parsed")

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
