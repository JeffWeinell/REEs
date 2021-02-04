#' Split character string to matrix
#'
#' Split a string into a matrix containing the splitted strings
#'
#' @param x Input character string to be split
#' @param split Character delimeter to split by (default is to split at spaces)
#' @param byrow Whether values in the output matrix should be filled by rows (default TRUE) or columns
#' @return Output matrix
#' @export
mat.strsplit <- function(x,split=" ",byrow=T){
	
	split.pattern <- split
	step1         <- strsplit(x,split=split.pattern) ### a list of vectors
	
	if(byrow==T){
		result <- do.call(rbind,step1)
	}
	if(byrow==F){
		result <- do.call(cbind,step1)
	}
	result
}
