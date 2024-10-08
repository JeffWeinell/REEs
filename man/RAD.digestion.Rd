% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddRADlike.R
\name{RAD.digestion}
\alias{RAD.digestion}
\title{SbfI.Seq                              <- REEs::datasets("restriction.enzymes")[datasets("restriction.enzymes")\link{,"RE.name"}=="SbfI","recognition.sequence"]
EcoR1.Seq                             <- REEs::datasets("restriction.enzymes")[datasets("restriction.enzymes")\link{,"RE.name"}=="EcoRI","recognition.sequence"]
Thermophis.baileyi_genome.url         <- REEs::datasets(1)[which(datasets(1)\link{,"species"}=="Thermophis baileyi"),"genome.url"]
proposed.ddRAD.loci.coordinates       <- REEs::proposeLoci.ddRADlike(input.seqs=Thermophis.baileyi_genome.url,output.dir="/ddRAD-like/",recognitionSeqs=c(SbfI.Seq,EcoR1.Seq),lim.lengths=c(900,1000),save.tables=T)
In Silico ddRAD version2}
\usage{
RAD.digestion(
  input.seqs,
  output.path = NULL,
  recognitionSeq,
  cut.after,
  include.antisense = T,
  include.uncut.contigs = T,
  fragment.size.filter = NULL
)
}
\arguments{
\item{input.seqs}{A DNAStringSet object or a character string indicating a local filepath or URL path to the input sequences in fasta format.}

\item{output.path}{Character string indicating where to save the output. If NULL (the default), the output is returned as an object but not saved to a file.}

\item{recognitionSeq}{A character vector containing the recognition sequence of each restriction enzyme.}

\item{cut.after}{Numerical vector of length ≥ 1, indicating the position of each cut site within each recognition sequence. A value of zero indicates that the cut site is immediately before the first base of the recognition sequence.}

\item{include.antisense}{Should antisense fragments also be include. Default TRUE. Probably should keep as TRUE for restriction enzymes that produce sticky-ends.}

\item{include.uncut.contigs}{Should contigs that did not contain any recognition sites be included in the output. Default is TRUE.}

\item{fragment.size.filter}{Either NULL (default) or a numerical vector of length 2 that defines the minimum and maximum length fragments to keep. If NULL, sequences are not filtered by length.}

\item{save.tables}{Should the tables holding the locations of each recognition sequence on each contig be saved? Default is TRUE but these tables can be big.}
}
\value{
A DNAStringSet object containing digested DNA sequences.
}
\description{
Perform in silico ddRAD on a reference genome.
This should be closer to actual ddRAD
}
