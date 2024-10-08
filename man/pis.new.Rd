% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pis.new.R
\name{pis.new}
\alias{pis.new}
\title{Parsimony Informative Sites}
\usage{
pis.new(x, abs = TRUE, use.ambiguities = FALSE, as.percent = F)
}
\arguments{
\item{x}{input alignment of DNA sequences stored as either a DNAStringSet, DNAMultipleAlignment, or DNAbin. This is an improved version of the pis function from the phyloch package.}

\item{abs}{Should the output of this function be the absolute number of parsimony informative sites? Default is TRUE. If FALSE, the argument as.percent determines the output type.}

\item{use.ambiguities}{Should ambiguous bases be considered as non-missing data? Default is false (Don't change this).}

\item{as.percent}{When abs = FALSE, should the output of this function be a percent? Default is FALSE, in which case the output is the fraction of sites parsimony informative.}
}
\value{
If abs = TRUE, this function returns the number of parsimony informative sites in the alignment.
If abs=FALSE and as.percent = TRUE, this function returns the percent of sites that are parsimony informative.
If abs=FALSE and as.percent = FALSE, this function returns the fraction of sites that are parsimony informative.
}
\description{
Calculates either the number of parsimony informative sites or the fraction or percent of sites that parsimony inforative.
}
