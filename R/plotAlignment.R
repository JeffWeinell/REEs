#' Plot Alignment
#' 
#' Plots the DNA alignment in the plotting window of R. Sequences are shown as lines where data is non-missing.
#' 
#' @param alignment Input DNA alignment of class DNAStringSet.
#' @param title Character string to use for the title. Default is no title.
#' @param colors Character vector with colors to use for 17 recognized characters in DNA alignments. If a single value is supplied, all non-gap characters are plotted with the sample supplied color.
#' Default colors "standard" = c(A="red", C="blue", G="yellow", T="green", M="darkgray", R="darkgray", W="darkgray", S="darkgray", Y="darkgray", K="darkgray", V="darkgray", H="darkgray", D="darkgray", B="darkgray", N="darkgray", '-'="white", '?'="white")
#' Note: "-" and "?" characters cannot be changed from background plot color.
#' @return An object containing the plot that is drawn in the plotting window of R.
#' @export
plotAlignment <- function(alignment,title="",colors="standard"){
	nsamples  <- length(alignment)
	width.al  <- width(alignment[1])
	#xvalsA    <- seq(from=1,to=(width.al+100),by=100)
	xvalsA    <- seq(from=1,to=(width.al),by=1)
	xvalsB    <- rep(xvalsA,nsamples)
	yvals     <- rep(c(1:nsamples),length(xvalsA))
	if(class(alignment)=="DNAStringSet"){
		CHARS=paste0(c(names(IUPAC_CODE_MAP),"-","\\?"),"+")
		CHARS_COLS=c(A="red",C="blue",G="yellow",T="green",M="darkgray", R="darkgray", W="darkgray", S="darkgray", Y="darkgray", K="darkgray", V="darkgray", H="darkgray", D="darkgray", B="darkgray", N="darkgray", '-'="white", '?'="white")
	}
	if(class(alignment)=="AAStringSet"){
		CHARS=paste0(c(names(AMINO_ACID_CODE),"-","\\?"),"+")
		CHARS_COLS <- c(D="red",E="red",C="yellow",M="yellow",U="yellow",K="blue",O="blue",R="blue",S="orange",T="orange",F="darkblue",Y="darkblue",N="cyan",Q="cyan",G="lightgray",L="green",V="green",I="green",A="lightgray",W="pink",H="paleblue",P="flesh",B="gray",J="gray",Z="gray",X="darkgray", '-'="white", '?'="white")
	}
	if(all(colors!="standard")){
		if(length(colors)==1 && is.null(names(colors))){
			CHARS_COLS[1:(length(CHARS_COLS)-2)] <- colors
		} else {
			CHARS_COLS[(names(CHARS_COLS) %in% names(colors))]  <- colors
			CHARS_COLS[!(names(CHARS_COLS) %in% names(colors))] <- "white"
		}
	}
	segments.list <- list(); length(segments.list) <- length(CHARS_COLS)
	for(i in 1:length(CHARS_COLS)){
		segments.x   <- stringr::str_locate_all(alignment,pattern=CHARS[i])
		nsegments    <- (lengths(segments.x)/2)
		sumnsegments <- sum(nsegments)
		segments.y   <- lapply(1:length(segments.x),function(x) {matrix(data=x,ncol=2,nrow=nsegments[x])})
		segments.xy  <- as.data.frame(cbind(do.call(rbind,segments.x), do.call(rbind,segments.y)))
		colnames(segments.xy) <- c("x0","x1","y0","y1")
		segments.xy[,"char"] <- rep(names(CHARS_COLS[i]),sumnsegments)
		segments.xy[,"char_color"] <- rep(unname(CHARS_COLS[i]),sumnsegments)
		segments.list[[i]] <- segments.xy
	}
	segments.df   <- do.call(rbind,segments.list)
	segments.mat  <- as.matrix(segments.df[,1:4])
	segments.cols <- as.matrix(segments.df[,6],drop=T)
	mat.plot      <- segments.mat[(segments.df[,"char"] %in% names(CHARS_COLS)),]
	mat.plot[,1]  <- (mat.plot[,1])-0.5
	mat.plot[,2]  <- (mat.plot[,2])+0.5
	mat.plot.cols <- segments.cols[(segments.df[,"char"] %in% names(CHARS_COLS))]
	plot(range(xvalsB),range(yvals),col="white",main=title,xlab="position",ylab="sample",ylim = rev(range(yvals)))
	rug(x = 1:nsamples, ticksize = -0.01, side = 2,quiet=T)
	segments(x0=mat.plot[,"x0"],x1=mat.plot[,"x1"],y0=mat.plot[,"y0"],y1=mat.plot[,"y1"],col=mat.plot.cols,lwd=2,lend="butt")
}

