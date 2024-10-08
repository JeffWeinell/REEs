% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRequired.R
\name{package.load.order}
\alias{package.load.order}
\title{State Required Packages}
\usage{
package.load.order(pkgs, search.for.more = T, order.by.repos = F)
}
\arguments{
\item{pkgs}{Character vector of package names.}

\item{search.for.more}{Should the output include dependencies of dependencies, etc.? Default TRUE.}

\item{order.by.repos}{Not yet implemented.}
}
\value{
A list of character vectors with the package names. Index of list indicates relative order of groups of packages that should be installed or loaded.
}
\description{
Get a list of all packages required (imports and depends) for a vector of packages.
Package versions are ignored.
}
