#' Locate Xth Occurrence in String
#' 
#' Character position of the Xth occurrence of a pattern. If the pattern is >1 character long, the function returns the position of the last character of the Xth occurrence.
#' 
#' @param strings One or more character strings to search within
#' @param pattern Character string to search for
#' @param X Integer indicating which occurrence to find
#' @return Integer vector indicating position of last occurrence of pattern in strings.
#' @export
str_locate_X <- function(strings,pattern,X){
 	all.locs <- stringr::str_locate_all(strings,pattern)
 	X.locs <- matrix(ncol=2,nrow=length(strings))
 	for(i in 1:length(strings)){
 		X.locs[i,1]      <- all.locs[[i]][X,1]
 		X.locs[i,2]      <- all.locs[[i]][X,2]
 	} 	
 	if(all(X.locs[,1]==X.locs[,2])){
 		X.locs <- c(X.locs[,1])
 	}
 	X.locs
}
