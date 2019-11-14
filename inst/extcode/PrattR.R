#' Automatically run Pratt and parse the results
#'
#' @param sequences A vector of sequences.
#' @param api The path to pratt.py API.
#' @param minPerc the minimum percentage of the input sequences that should match a pattern (C%).
#' @param maxPatternLength Maximum pattern length (PL parameter) allows you to set the maximum length of a pattern.
#' @param maxNumPatternSymbols Maximum number of pattern symbols (PN parameter).
#' @param maxNumWildcard Maximum length of a widecard (x).
#' @param maxNumFlexSpaces Maximum number of flexible space (matching a variable number of arbitrary sequence symbols).
#' @param maxFlexibility The maximum flexibility of a flexible wildcard (matching a variable number
#' of arbitrary sequence symbols) (FL parameter). For instance x(2,4) and x(10,12) has flexibility 2.
#' @param maxFlexProduct Upper limit on the product of a flexibilities for a pattern.
#' @param patternRefinement Pattern Refinement (R parameter).
#' @param patternFormat PROSITE Pattern Format (OP parameter). When switched on, patterns
#' will be output in PROSITE style (for instance C-x(2,4)-[DE]).
#' @param maxNumPatterns Maximum number of patterns (ON parameter) between 1 and 100.
#' @param printVertically Print vertically (MV parameter). if set, the output is printed vertically instead
#' of horizontally, vertical output can be better for large sequence sets.
#'
#' @export
#' @importFrom data.table fread
#'
PrattR <- function(sequences, api = "/Users/Wubing/Jobs/Project/UPS/_Code/pratt.py",
                   minPerc = 20, maxPatternLength = 15,
                   maxNumPatternSymbols = maxPatternLength,
                   maxNumWildcard = 9, maxNumFlexSpaces = 7,
                   maxFlexibility = 2, maxFlexProduct = 50,
                   patternRefinement = TRUE, patternFormat = "PROSITE",
                   maxNumPatterns = 100, printVertically = TRUE){
  requireNamespace("data.table")
  if(is.null(names(sequences)))
    names(sequences) = paste0("inputseq_", 1:length(sequences))
  if(length(sequences)>100) sequences = sequences[1:100]
  tmp =  rep(sequences, each = 2)
  tmp[seq(1, length(tmp), 2)] = paste0(">", names(sequences))
  write.table(tmp, "tmp.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
  jobid = system(paste0("/Users/Wubing/miniconda2/bin/python ", api,
                        " --asyncjob --email zwbing_dyup@tongji.edu.cn",
                        " --minPerc ", minPerc,
                        " --maxPatternLength ", maxPatternLength,
                        " --maxNumWildcard ", maxNumWildcard,
                        " --maxNumPatternSymbols ", maxNumPatternSymbols,
                        " --maxNumFlexSpaces ", maxNumFlexSpaces,
                        " --maxFlexibility ", maxFlexibility,
                        " --maxFlexProduct ", maxFlexProduct,
                        ifelse(patternFormat=="PROSITE", " --patternFormat", ""),
                        ifelse(printVertically, " --printVertically", ""),
                        " --maxNumPatterns ", maxNumPatterns,
                        ifelse(patternRefinement, " --patternRefinement", ""),
                        " --sequence tmp.fasta"), intern = TRUE)
  file.remove("tmp.fasta")

  ## Read output file
  host = "https://www.ebi.ac.uk/Tools/services/rest/pratt/result/"
  while (TRUE) {
    status = system(paste0("/Users/Wubing/miniconda2/bin/python ", api, " --status --jobid ", jobid[1]), intern = TRUE)
    if(status[2]=="FINISHED"){
      out = data.table::fread(paste0(host, jobid[1], "/out"), fill = TRUE)
      break
    }
  }

  ## Prioritize the output mode
  idx1 = which(out[,1]=="Best Patterns before refinement:")
  idx2 = which(out[,1]=="Best Patterns (after refinement phase):")
  idx3 = which(out[,1]=="Best patterns with alignments:")
  refine_before = out[(idx1+2):(idx2-3), ]
  tmp = lapply(strsplit(unlist(refine_before), split = " "), function(x) x[x!=""])
  refine_before = t(as.data.frame(lapply(tmp, function(x) x[x!=""])))
  refine_before = refine_before[, -1]
  colnames(refine_before) = c("Fitness", "Hits", "Seqs", "Pattern")
  refine_before[,2] = gsub("\\(", "", refine_before[,2])
  refine_before[,3] = gsub("\\)", "", refine_before[,3])
  rownames(refine_before) = 1:nrow(refine_before)
  refine_before = as.data.frame(refine_before, stringsAsFactors = FALSE)
  refine_before$Fitness = as.numeric(refine_before$Fitness)
  refine_before$Hits = as.integer(refine_before$Hits)
  refine_before$Seqs = as.integer(refine_before$Seqs)

  refine_after = out[(idx2+2):(idx3-3), ]
  tmp = strsplit(gsub("^[A-z] *", "", unlist(refine_after)), split = " ")
  tmp = lapply(tmp, function(x){ x[x!=""]})
  refine_after = t(as.data.frame(tmp))
  refine_after = refine_after[, -1]
  colnames(refine_after) = c("Fitness", "Hits", "Seqs", "Pattern")
  refine_after[,2] = gsub("\\(", "", refine_after[,2])
  refine_after[,3] = gsub("\\)", "", refine_after[,3])
  rownames(refine_after) = 1:nrow(refine_after)
  refine_after = as.data.frame(refine_after, stringsAsFactors = FALSE)
  refine_after$Fitness = as.numeric(refine_after$Fitness)
  refine_after$Hits = as.integer(refine_after$Hits)
  refine_after$Seqs = as.integer(refine_after$Seqs)

  res = list(outsite = paste0(host, jobid[1], "/out"),
             rawoutput = out,
             refine_before = refine_before,
             refine_after = refine_after)
  return(res)
}
