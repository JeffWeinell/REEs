#' Plot Alignment
#' 
#' Plots the DNA alignment in the plotting window of R. Sequences are shown as lines where data is non-missing.
#' 
#' @param alignment Input DNA alignment of class DNAStringSet.
#' @return An object containing the plot that is drawn in the plotting window of R.
#' @export
plotAlignment <- function(alignment){
	nsamples  <- length(alignment)
	width.al  <- width(alignment[1])
	xvalsA    <- seq(from=1,to=(width.al+100),by=100)
	xvalsB    <- rep(xvalsA,nsamples)
	yvals     <- rep(c(1:nsamples),length(xvalsA))
	plot(xvalsB,yvals,col="white",xlab="position",ylab="sample",ylim = rev(range(yvals)))
	for(i in 1:length(alignment)){
		gapLocation   <- unique(unlist(stringr::str_locate_all(alignment[i],pattern="-")))
		if(length(gapLocation)==0){
			noGapLocation <- c(1:width.al)
		} else {
			noGapLocation <- c(1:width.al)[-gapLocation]
		}
		points(noGapLocation,rep(i,length(noGapLocation)),pch=15,cex=0.5)
	}
}
