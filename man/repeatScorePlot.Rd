% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexity.R
\name{repeatScorePlot}
\alias{repeatScorePlot}
\title{Plot sliding window DNA repeat score

Plots sequence repeat score for a specified repeat pattern along an input DNA sequence along sliding windows and for a random sequence of the same length.

Black line shows repeat score of the input DNA sequence for a sliding window with windowsize=10 and windowshift=5
Red line shows repeat score of the indput DNA sequence for a sliding window with windowsize=50 and windowshift=50
Blue line shows repeat score of the randomly generated DNA sequence using a sliding window with windowsize=50 and windowshift=50
Note: I may change this to a ggplot and provide options to control sliding windows}
\usage{
repeatScorePlot(dna, rp = "1.1", wsize = 1)
}
\arguments{
\item{dna}{Character string, or a DNAString object, or a length 1 DNAStringSet object as defined by the Biostrings package.}

\item{rp}{Text string defining that defines how complexity is defined. \link{AT PRESENT ONLY PATTERNS COMPOSED OF "1" and "." CAN BE USED}. Default is "1.1", which specifies to return the fraction of bases in which base\link{i} == base\link{i+2}.}

\item{wsize}{Number that controls windowsizes; default window sizes of the larger windows (indicated by red and blue lines) are multiplied by this number.}
}
\description{
Plot sliding window DNA repeat score

Plots sequence repeat score for a specified repeat pattern along an input DNA sequence along sliding windows and for a random sequence of the same length.

Black line shows repeat score of the input DNA sequence for a sliding window with windowsize=10 and windowshift=5
Red line shows repeat score of the indput DNA sequence for a sliding window with windowsize=50 and windowshift=50
Blue line shows repeat score of the randomly generated DNA sequence using a sliding window with windowsize=50 and windowshift=50
Note: I may change this to a ggplot and provide options to control sliding windows
}
