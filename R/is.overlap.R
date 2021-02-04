#' @title Test if Two Ranges Overlap
#'
#' This function takes as input two integer vectors. The min and max values of each vector represent the lower and upper limits of each range, respectively. Returns TRUE if input ranges overlap, otherwise FALSE.
#'
#' @param r1 Integer vector with min and max values representing limits of first range
#' @param r2 Integer vector with min and max values representing limits of second range
#' @return Logical indicating whether or not r1 and r2 overlap
#' @export is.overlap
is.overlap <- function(r1,r2){
	 min(r1) <= max(r2) && min(r2) <= max(r1)
}


#' @title Test if subset
#' 
#' Test if a range of numbers is a subset of another range of numbers.
#' 
#' @param r1 Integer vector with min and max values representing limits of first range
#' @param r2 Integer vector with min and max values representing limits of second range
#' @return Logical indicating whether or not r1 is a non-equivalent subset of r2.
#' @export is.subset
is.subset <- function(r1,r2){
	test1  <- (min(r1) >= min(r2) && max(r1) <= max(r2))
	test2  <- any(c(min(r1),max(r1)) != c(min(r2),max(r2)))
	result <- all(test1,test2)
	result
}