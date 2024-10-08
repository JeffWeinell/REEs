% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trim.alignment.R
\name{trim.alignment}
\alias{trim.alignment}
\title{Trim DNA alignment}
\usage{
trim.alignment(alignment, threshold = 1)
}
\arguments{
\item{alignment}{DNA alignment held in an object of class DNAStringSet, DNAMultipleAlignment, or DNAbin. These object classes are defined in the Biostrings and ape packages.}

\item{threshold}{Number of individuals with non-missing data that must be present to keep the alignment ends. Default is 1.}
}
\value{
Trimmed alignment object having the same class as the input alignment.
}
\description{
Trims the first and last columns of a DNA alignment until those columns contain at least the defined threshold number of individuals with non-missing data.
}
