% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codonsDNAStringSet.R
\name{codonsDNAStringSet}
\alias{codonsDNAStringSet}
\title{Get Sites at Codon Position}
\usage{
codonsDNAStringSet(
  x,
  which.codon.positions,
  predict.reading.frame = F,
  frameAtFirstchar = NULL
)
}
\arguments{
\item{x}{Input alignment DNAStringSet object.}

\item{which.codon.positions}{An integer in {1,2,3}, or an integer vector containing the values 1, 2, and/or 3.}

\item{predict.reading.frame}{Whether or not to try to guess the reading frame, based on the relative number of snps in each frame. Default is FALSE.
The third frame is assumed to be the frame with the most snps. This method works best for longer, faster-evolving loci, and when the alignment has many individuals, and/or moderate to high genetic distances among individuals.}

\item{frameAtFirstchar}{The frame (i.e., codon position) of the first character of the sequence (if known). If unknown, leave as NULL, in which case the first base treated as if it is the first codon position.}
}
\value{
DNAStringSet object holding only the sites at codon positions indicated by the argument which.codon.positions
}
\description{
Returns a DNA alignment of only the first, second, or third codon position sites (or any two of these three).
}
