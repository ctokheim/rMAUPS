% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsScore.R
\docType{methods}
\name{gsScore}
\alias{gsScore}
\title{Calculate score across genes and samples}
\usage{
gsScore(dat, gset, fun = "PC")
}
\arguments{
\item{dat}{normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples}

\item{gset}{Gene sets.}

\item{fun}{("PC", default), Pearson, ssGSEA or mean (other value). fisher, stouffer}
}
\value{
numeric vector or matrix.
}
\description{
This wrapper function combines filtering out genes with low reads in a number of
samples (recommended for limma:voom) with normalization
}
\author{
Wubing Zhang
}
