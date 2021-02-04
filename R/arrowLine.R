#' Plot Arrow
#'
#' Draw a line with arrows in it; e.g. -->-->-->
#' This was from a comment by Stack Overflow user Sacha Epskamp
#'
#' @param x0 x start position
#' @param y0 y start position
#' @param x1 x end position
#' @param y1 y end position
#' @param nArrow arrow length
#' @return 
#' @export 
arrowLine <- function(x0,y0,x1,y1=y0,nArrow=1,...){
	lines(c(x0,x1),c(y0,y1),...)
	Ax=seq(x0,x1,length=nArrow+1)
	Ay=seq(y0,y1,length=nArrow+1)
	for (i in 1:nArrow){
		arrows(Ax[i],Ay[i],Ax[i+1],Ay[i+1],...)
	}
}
