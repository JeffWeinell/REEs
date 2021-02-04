#' Split Columns by Delimiter
#'
#' Split a matrix column into multiple columns by a delimiter
#'
#' @param x Input matrix and column to be split
#' @param split Character to delimit columns by
#' @param ncol number of columns in final matrix (default is to determine automatically)
#' @param byrow Logical indicating if values in the result matrix should be filled by rows and then columns
#' @return Output matrix
#' @export
column.split <- function(x,split=" ",ncol="auto",byrow=T){
	if(is.vector(x,mode="character")){
		x <- matrix(x,ncol=1)
	}
	result.temp <- apply(X=x,MARGIN=1,FUN=mat.strsplit)
	result      <- do.call(rbind, result.temp)
	result
}
