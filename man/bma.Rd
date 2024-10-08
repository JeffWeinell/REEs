% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bma.R
\name{bma}
\alias{bma}
\title{IUPAC nucleotide symbols

This function returns the IUPAC symbol for a nucleotide sequence, for instance c("c", "c", "g") is coded by "s",
which is useful for describing ambiguous or consensus nucleotide sequences.
This function is identical to the function bma in package seqinr.}
\usage{
bma(nucl, warn.non.IUPAC = TRUE, type = c("DNA", "RNA"))
}
\arguments{
\item{nucl}{A nucleotide sequence as a vector of single characters (e.g., nucl=c("c", "c", "g")).}

\item{warn.non.IUPAC}{If TRUE (default) warns when no IUPAC symbol is possible.}

\item{type}{Type of nucleotides: either "DNA" or "RNA"}
}
\value{
A single IUPAC symbol in lower case, or NA when this is not possible.
}
\description{
IUPAC nucleotide symbols

This function returns the IUPAC symbol for a nucleotide sequence, for instance c("c", "c", "g") is coded by "s",
which is useful for describing ambiguous or consensus nucleotide sequences.
This function is identical to the function bma in package seqinr.
}
