#' Intersect of all sets
#'
#' Get the intersect set of two or more sets
#'
#' @param x Input list of vectors
#' @return Output vector of elements shared among x
#' @export
intersect.all <- function(x){
	result <- Reduce(intersect,x)
	result
}
