#' Matrix to List of Strings
#' 
#' Collapses rows of a character matrix into character strings
#' 
#' @param y Input character matrix
#' @return  Returns a list character strings, each string formed by the row elements of the input matrix
#' @export 
string.list <- function(y){
	paste.collapse <- function(z){
		z <- paste(z,collapse="")
		z
	}
	if(class(y)[1]!="matrix"){
		y <- as.matrix(y)
	}
	y <- apply(X=y,MARGIN=1,FUN=paste.collapse)
	y
}
