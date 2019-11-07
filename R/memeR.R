#' Run MEME using R
#'
#' @param sequences a list of sequences.
#' @param background a list of control sequences.
#' @param outdir name of directory for output files.
#' @param objfun objective function (default: classic).
#' Classic mode: You provide one set of sequences and MEME discovers motifs
#' enriched in this set. Enrichment is measured relative to a (higher order)
#' random model based on frequencies of the letters in your sequences, or relative
#' to the frequencies given in a "Custom background model" that you may provide.
#' Discriminative mode: You provide two sets of sequences and MEME discovers motifs
#' that are enriched in the first (primary) set relative to the second (control) set.
#' In Discriminative mode, we first calculate a position-specific prior from the two
#' sets of sequences. MEME then searches the first set of sequences for motifs using
#' the position-specific prior to inform the search. This approach is based on the
#' simple discriminative prior "D" described in Section 3.5 of Narlikar et al.
#' We modify their approach to search for the "best" initial motif width, and to handle
#' protein sequences using spaced triples. Refer to the psp-gen documentation and to
#' our paper for more details.
#' Differential Enrichment mode: You provide two sets of sequences and MEME discovers
#' motifs that are enriched in the first (primary) set relative to the second (control)
#' set. In Differential Enrichment mode, MEME optimizes an objective function based on
#' the hypergeometric distribution to determine the relative enrichment of sites in the
#' primary sequences compared to the control sequences.
#' @param test statistical test type (default: mhg).
#' @param shuf preserve frequencies of k-mers of size <kmer> when shuffling.
#' @param seed random seed for shuffling and sampling.
#' @param hsfrac fraction of primary sequences in holdout set.
#' @param cefrac fraction sequence length for CE region.
#' @param searchsize maximum portion of primary dataset to use for motif search (in characters).
#' @param alph sequences use DNA/RNA/protein alphabet.
#' @param revcomp allow sites on + or - DNA strands.
#' @param mod distribution of motifs.
#' @param nmotifs maximum number of motifs to find.
#' @param evt stop if motif E-value greater than <evt>.
#' @param time quit before <t> CPU seconds consumed.
#' @param minsites minimum number of sites for each motif.
#' @param maxsites maximum number of sites for each motif.
#' @param minw minimum motif width.
#' @param maxw maximum motif width.
#' @param opencost gap opening cost for multiple alignments.
#' @param extendcost gap opening cost for multiple alignments.
#' @param noendgaps do not count end gaps in multiple alignments.
#' @param maxiter maximum EM iterations to run.
#' @param prior type of prior to use.
#'
#' @export
#' @importFrom data.table fread
#'
memeR <- function(sequences, background = NULL, outdir = "./",
                  objfun = c("classic", "de", "se", "cd", "ce")[1],
                  test = c("mhg", "mbn", "mrs")[1],
                  shuf = 10, seed = 50, hsfrac = 0.5, cefrac = 0.25,
                  searchsize = 0.5, alph = "protein", revcomp = FALSE,
                  alphfile = "/Users/Wubing/Jobs/Project/UPS/_Data/AA_letters",
                  mod = c("oops", "zoops", "anr")[1],
                  nmotifs = 5, evt = NULL, time = 6000,
                  minsites = 0.2*length(sequences), maxsites = 600,
                  minw = 6, maxw = 15,
                  opencost = NULL, extendcost = NULL, noendgaps = TRUE, maxiter = 30,
                  prior = c("dirichlet", "dmix", "mega", "megap", "addone")[1],
                  api = "/Users/Wubing/Applications/meme/bin/meme"){
  requireNamespace("data.table")
  if(is.null(names(sequences)))
    names(sequences) = paste0("primary_", 1:length(sequences))
  tmp =  rep(sequences, each = 2)
  tmp[seq(1, length(tmp), 2)] = paste0(">", names(sequences))
  write.table(tmp, "primary.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(!is.null(background)){
    if(is.null(names(background)))
      names(background) = paste0("background_", 1:length(background))
    tmp = rep(background, each = 2)
    tmp[seq(1, length(tmp), 2)] = paste0(">", names(background))
    write.table(tmp, "background.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  joblog = system(paste0(api, " primary.fasta -o ", outdir,
                        " -objfun ", objfun,
                        " -shuf ", shuf, " -seed ", seed,
                        ifelse(objfun=="classic", "", paste0(
                          " -test ", test,
                          ifelse(is.null(background), "", " -neg background.fasta"),
                          " -hsfrac ", hsfrac,
                          ifelse(objfun=="ce", paste0(" -cefrac ", cefrac), ""),
                          " -searchsize ", searchsize)),
                        " -", alph, " -alph ", alphfile, ifelse(revcomp, " -revcomp", ""),
                        " -mod ", mod, " -nmotifs ", nmotifs,
                        ifelse(is.null(evt), "", paste0(" -evt ", evt)),
                        " -time ", time, " -minsites ", minsites,
                        " -maxsites ", maxsites, " -minw ", minw, " -maxw ", maxw,
                        ifelse(is.null(opencost), "", paste0(" -wg ", opencost)),
                        ifelse(is.null(extendcost), "", paste0(" -ws ", extendcost)),
                        ifelse(noendgaps, " -noendgaps", ""),
                        " -maxiter ", maxiter, " -prior ", prior), intern = TRUE)
  # suppressWarnings(file.remove(c("primary.fasta", "background.fasta")))
  meme_res = parseMEME(outdir = outdir)
  return(meme_res)
}


parseMEME <- function(outdir = "testMEME/"){
  meme_raw = data.table::fread(paste0(outdir, "meme.txt"),
                               stringsAsFactors = FALSE,
                               fill = TRUE, skip = 25)
  meme_raw = as.data.frame(meme_raw, stringsAsFactors = FALSE)
  ALPHABET = meme_raw$V2[meme_raw$V1=="ALPHABET="]
  idx1 = which(meme_raw$V1=="Sequence")[1]+2
  idx2 = which(grepl("^\\*", meme_raw$V1))
  idx2 = idx2[idx2>idx1][1] - 1
  sequences = unlist(meme_raw[idx1:idx2, seq(1, ncol(meme_raw),3)])
  sequences = sequences[!(sequences=="" | is.na(sequences))]

  idx1 = which(grepl("Motif", meme_raw$V1) & grepl("^MEME", meme_raw$V3) & meme_raw$V4=="Description")
  motifsummary = meme_raw[idx1-3, ]
  motifsummary = motifsummary[, c(2,5,8,11,14)]
  colnames(motifsummary) = c("Motif", "Width", "Sites", "Llr", "E-value")
  rownames(motifsummary) = motifsummary$Motif
  MotifInSeq <- lapply(motifsummary$Motif, function(m, meme_raw){
    idx1 = which(grepl("Motif", meme_raw$V1) & grepl(m, meme_raw$V2) & grepl("sorted", meme_raw$V5))
    tmp = meme_raw[(idx1+4):(length(sequences)+idx1), 1:6]
    colnames(tmp) = c("Sequence", "Start", "P-value", "Before", "Site", "After")
    tmp = cbind(Motif = m, tmp)
  }, meme_raw)
  names(MotifInSeq) = motifsummary$Motif
  res = list(Sequences = sequences, MotifSummary = motifsummary, MotifInSeq = MotifInSeq)
  return(res)
}


