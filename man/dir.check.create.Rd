% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dir.check.create.R
\name{dir.check.create}
\alias{dir.check.create}
\title{Create directory unless it already exists}
\usage{
dir.check.create(directory)
}
\arguments{
\item{directory}{name or path to the directory to check and create}
}
\value{
A new directory unless it already existed
}
\description{
Checks if a directory exists, and if not, creates it. Parent directories are also created if they do not already exist.
}
