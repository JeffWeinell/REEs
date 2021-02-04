#' Locate Last in Strings
#' 
#' Character position of the last occurrence of pattern if the pattern is >1 character long, the function returns the position of the last character of the last instance of the pattern
#' 
#' @param strings One or more character strings to search within
#' @param pattern Character string to search for
#' @return Integer vector indicating position of last occurrence of pattern in strings.
#' @export
str_locate_last <- function(strings,pattern){
 	all.locs  <- stringr::str_locate_all(strings,pattern)
 	last.locs <- list()
 	length(last.locs) <- length(strings)
 	for(i in 1:length(strings)){
 		last.locs[i] <- max(all.locs[[i]][,2])
 	}
 	last.locs <- unlist(last.locs)
 	last.locs
}
