% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.mafft.R
\name{run.mafft}
\alias{run.mafft}
\title{Run MAFFT}
\usage{
run.mafft(
  mafft.path = NULL,
  unaligned.contigs,
  add.contigs = NULL,
  algorithm = "localpair",
  rev.dir = T,
  save.name = NULL,
  threads = 6,
  delete.files = T
)
}
\arguments{
\item{mafft.path}{Path to MAFFT executable. Default NULL, in which case the mafft executable is searched for on your path and environement.}

\item{unaligned.contigs}{Input DNA sequences in an object of class DNAStringSet.}

\item{add.contigs}{Default is NULL. This argument not currently implemented in this function.}

\item{algorithm}{Alignment algorithm to use. Default is "localpair". Alternatives are "genafpair" and "globalpair". See MAFFT documentation for more details: https://mafft.cbrc.jp/alignment/software/manual/manual.html.}

\item{rev.dir}{Should MAFFT try to align reverse complement of sequences? Default is TRUE.}

\item{save.name}{Filename for DNA alignment. This file is temporary unless the argument delete.files is changed to FALSE.}

\item{threads}{How many threads should be used (default is 6)}

\item{delete.files}{If files should be deleted (default is T)}
}
\value{
DNAStringSet object holding the aligned DNA sequences
}
\description{
R wrapper for running MAFFT alignment software.
Need to check if there is any reason to use this function instead of the REEs::mafft
}
