#' Replace multiple patterns
#' 
#' Replace a set of characters with a set of different characters
#'
#' @param x A character vector of patterns to be replaced
#' @param y A character vector of replacement patterns
#' @param z Input vector
#' @return output with x values replaced by y values
#' @export 
mgsub <- function(x,y,z){
	for(i in 1:length(x)){
		z <- gsub(pattern=x[i],replacement=y[i],z)
	}
	z
}
