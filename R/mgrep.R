#' Find patterns in subjects
#' 
#' Find location of each value in a query vector within a subject vector
#'
#' @param query A character vector of patterns to find
#' @param subject A character vector of subjects to search within
#' @return A list numeric vectors. The ith vector indicates which elements of the subject vector contains the ith query.
#' @export 
mgrep <- function(query,subject){
	z <- NULL
	for(i in 1:length(query)){
		z[i] <- grep(pattern=query[i],x=subject)
	}
	z
}
