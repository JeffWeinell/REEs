#' Round in specific way
#'
#' Round up or down to the nearest supplied base value
#'
#' @param x input number to be rounded
#' @param base what number to round by (e.g., base = 5 if "round to the nearest 5")
#' @param direction "nearest" (default), alternatively "up" or "down"
#' @return rounded number
#' @export
mround <- function(x,base=1,direction="nearest"){
	res <- base*round(x/base)
	if(direction=="up" & res < x){
		res <- res+base
	}
	if(direction=="down" & res > x){
		res <- res-base
	}
	res
}
