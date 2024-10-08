% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexity.R
\name{complexity}
\alias{complexity}
\title{Sequence Complexity

This function calculates sequence complexity as the fraction of bases that differ from the next base (base\link{i} != base\link{i+1}). This is the same definition used in the program fastp: https://github.com/OpenGene/fastp#low-complexity-filter}
\usage{
complexity(x, w = 1)
}
\arguments{
\item{x}{Character string, or a DNAString object, or a length 1 DNAStringSet object as defined by the Biostrings package.}

\item{w}{The number of subsequent bases to compare to the ith base. Default 1.}
}
\value{
Number between 0 and 1 indicating the sequence complexity score. A value of 1 indicates that every base differs from the preceding w bases, where as a value of 0 indicates that all bases are the same.
}
\description{
Sequence Complexity

This function calculates sequence complexity as the fraction of bases that differ from the next base (base\link{i} != base\link{i+1}). This is the same definition used in the program fastp: https://github.com/OpenGene/fastp#low-complexity-filter
}
