#' Rolling Summation
#'
#' Returns the rolling summation vector of input numerical vector
#' 
#' @param x Input numerical vector
#' @return Output is a numerical vector that is the rolling summation of the input vector x
#' @export
rollSum <- function(x){
	rollingSum <- list()
	length(rollingSum) <- length(x)
	for(i in 1:length(x)){
		rollingSum[i] <- sum(x[1:i])
	}
	rollingSum <- unlist(rollingSum)
	rollingSum
}
