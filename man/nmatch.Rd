% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmatch.R
\name{nmatch}
\alias{nmatch}
\title{Find First N Matches}
\usage{
nmatch(x, y, n)
}
\arguments{
\item{x}{Vector of values to look for. x should not contain NA values.}

\item{y}{Vector to search within.}

\item{n}{Number of matches to return.}
}
\value{
when simplify = FALSE, this function returns a list of numeric vectors. The ith vector contains the locations of (at most) the first n matches of x\link{i} in y.
If simplify = TRUE, a single numeric vector is returned that corresponds to concatenating the list of vectors that would hav ebeen returned if simplify=F.
}
\description{
This function returns the indices of the first n matches of each value in a vector x present in a vector y.
}
