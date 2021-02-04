#' Count character occurences
#' 
#' Count the number of occurrences of a character in a character string
#' 
#' @param x character string to search within
#' @param pattern the character or character string to be counted in x
#' @return integer of the number of occurrence of pattern in x
#' @export
countChars <- function(x,pattern) {
	counter <- gsub(pattern,"",x)
	nchar(x)-nchar(counter)
}
