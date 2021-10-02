#' @title Sequence Complexity
#'
#' This function calculates sequence complexity as the fraction of bases that differ from the next base (base[i] != base[i+1]). This is the same definition used in the program fastp: https://github.com/OpenGene/fastp#low-complexity-filter
#' 
#' @param x Character string, or a DNAString object, or a length 1 DNAStringSet object as defined by the Biostrings package.
#' @param w The number of subsequent bases to compare to the ith base. Default 1.
#' @return Number between 0 and 1 indicating the sequence complexity score. A value of 1 indicates that every base differs from the preceding w bases, where as a value of 0 indicates that all bases are the same.
#' @export complexity
complexity <- function(x,w=1){
	if(class(x)=="DNAString"){
		x <- as.character(unlist(x))
	}
	if(length(x)>1){
		stop("Sequences must supplied as a DNAString object or a length 1 DNAStringSet object")
	}
	if(class(x)=="DNAStringSet"){
		x <- as.list(x[1])[[1]]
	}
	tempseq <- unlist(strsplit(as.character(unlist(x)),split=""))
	res     <- length(which(sapply(1:(length(tempseq)-1),function(i){length(unique(tempseq[i:(i+w)]))!=1})))/(length(tempseq)-1)
	res
}

#' @title Repeat Score
#'
#' This function calculates fraction of bases in a DNA sequence that initiate a specified pattern. At present only the default pattern (base[i] == base[i+2]) can be used.
#' 
#' @param x Character string, or a DNAString object, or a length 1 DNAStringSet object as defined by the Biostrings package.
#' @param rp Text string defining that defines how complexity is defined. [AT PRESENT ONLY PATTERNS COMPOSED OF "1" and "." CAN BE USED]. Default is "1.1", which specifies to return the fraction of bases in which base[i] == base[i+2].
#' @return Number between 0 and 1 indicating the repeat score for the repeat pattern indicated by rp. A value of 1 indicates that all possible reading frames and frame position having a width equal to the rp pattern width have the rp pattern. A value of zero indicates that the pattern is absent.
#' @export repeatScore
repeatScore <- function(x,rp="1.1"){
	if(class(x)=="DNAString"){
		x <- as.character(unlist(x))
	}
	if(length(x)>1){
		stop("Sequences must supplied as a DNAString object or a length 1 DNAStringSet object")
	}
	if(class(x)=="DNAStringSet"){
		x <- as.list(x[1])[[1]]
	}
	patternWidth <- nchar(rp)
	pattern      <- unlist(strsplit(rp,split=""))
	relpattern   <- which(pattern==pattern[1])
	tempseq      <- unlist(strsplit(as.character(unlist(x)),split=""))
	res          <- length(which(sapply(1:(length(tempseq)-(patternWidth-1)),function(i){length(unique(tempseq[relpattern+(i-1)]))==1})))/(length(tempseq)-(patternWidth-1))
	#if(rp=="1.1"){
	#	res     <- length(which(sapply(1:(length(tempseq)-2),function(i){tempseq[i]==tempseq[i+2]})))/(length(tempseq)-(patternWidth-1))
	#}
	res
}

#' @title Apply function as sliding window
#'
#' This function applies a function to a DNA sequence as a sliding window.
#' 
#' @param x Character string, or a DNAString object, or a length 1 DNAStringSet object as defined by the Biostrings package.
#' @param FUN The function to apply to the DNA sequence framed by each position of the window.
#' @param windowsize The size of the window in base pairs. Default value NULL, in which case window size is equal to the length of x.
#' @param windowshift The distance to move the window between steps. Default value NULL, in which case windowshift=windowsize.
#' @param startbase Number indicating the first base to include in the first frame position. Default 1.
#' @param remainder Logical indicating whether to apply the function to the last window if the window size is less than indicated by the windowsize argument. Default TRUE.
#' @param minwindowsize Number indicating the minimum allowable window size. Ignored if windowsize < minwindowsize.
#' @return A list or vector, with the ith entry equal to the value returned by applying FUN to the ith window of the input DNA sequence.
#' @export wapply
wapply <- function(x,FUN,windowsize=NULL,windowshift=NULL,startbase=1,remainder=T,minwindowsize=2){
	if(class(x)=="DNAString"){
		x <- as.character(unlist(x))
	}
	if(length(x)>1){
		stop("Sequences must supplied as a DNAString object or a length 1 DNAStringSet object")
	}
	if(class(x)=="DNAStringSet"){
		x <- as.list(x[1])[[1]]
	}
	tempseq       <- unlist(strsplit(as.character(unlist(x)),split=""))
	if(is.null(windowsize)){
		windowsize <- length(tempseq)
	}
	if(is.null(windowshift)){
		windowshift <- windowsize
	}
	window.starts <- seq(from=startbase,to=length(tempseq),by=windowshift)
	window.ends   <- window.starts+(windowsize-1)
	window.starts[!window.starts <= length(tempseq)] <- length(tempseq)
	window.ends[!window.ends <= length(tempseq)] <- length(tempseq)
	windows.mat   <- unique(cbind(window.starts,window.ends))
	if(!remainder){
		windows.mat <- windows.mat[((windows.mat[,2]-windows.mat[,1])+1)==windowsize,]
	}
	if(windowsize >= minwindowsize){
		windows.mat <- windows.mat[((windows.mat[,2]-windows.mat[,1])+1)>=minwindowsize,]
	}
	window.starts <- windows.mat[,1]
	window.ends   <- windows.mat[,2]
	funres <- list(); length(funres) <- length(window.starts)
	for(i in 1:length(window.starts)){
		windowseq <- DNAString(paste(tempseq[window.starts[i]:window.ends[i]],collapse=""))
		funres[[i]] <- FUN(windowseq)
	}
	classes <- lapply(funres, class)
	if(all(lengths(funres)==1) & length(unique(unlist(classes)))==1){
		funres <- unlist(funres)
	}
	funres
}

#' @title Plot sliding window DNA complexity
#'
#' Plots sequence complexity of an input DNA sequence along sliding windows and a random sequence of the same length.
#' Black line shows complexity of the input DNA sequence for a sliding window with windowsize=10 and windowshift=5
#' Red line shows complexity of the input DNA sequence for a sliding window with windowsize=200 and windowshift=100
#' Blue line shows complexity of the randomly generated DNA sequence using a sliding window with windowsize=200 and windowshift=100
#' Note: I may change this to a ggplot and provide options to control sliding windows
#' 
#' @param dna Character string, or a DNAString object, or a length 1 DNAStringSet object as defined by the Biostrings package.
#' @param wsize Number that controls windowsizes; default window sizes of the larger windows (indicated by red and blue lines) are multiplied by this number.
#' @return NULL
#' @export complexityPlot
complexityPlot <- function(dna,wsize=1){
	dna   <- DNAStringSet(dna)
	test1 <- wapply(x=dna, FUN=complexity, windowsize=10, windowshift=5)
	test2 <- wapply(x=dna, FUN=complexity, windowsize=200*wsize, windowshift=100*wsize)
	sim1  <- DNAStringSet(paste(sample(c("A","C","G","T"),size=width(dna),replace=T),collapse=""))
	test3 <- wapply(x=sim1, FUN=complexity, windowsize=200*wsize, windowshift=100*wsize)
	plot(seq(from=1,to=width(dna),length.out=length(test1)),test1,type="l",xlab="nucleotide position along sequence",ylab="DNA sequence complexity",sub="black and red = empirical; blue = random",ylim=c(0,1))
	lines(seq(from=1,to=width(dna),length.out=length(test2)),test2,col="red")
	lines(seq(from=1,to=width(sim1),length.out=length(test3)),test3,col="blue")
}


#' @title Plot sliding window DNA repeat score
#'
#' Plots sequence repeat score for a specified repeat pattern along an input DNA sequence along sliding windows and for a random sequence of the same length.
#' 
#' Black line shows repeat score of the input DNA sequence for a sliding window with windowsize=10 and windowshift=5
#' Red line shows repeat score of the indput DNA sequence for a sliding window with windowsize=50 and windowshift=50
#' Blue line shows repeat score of the randomly generated DNA sequence using a sliding window with windowsize=50 and windowshift=50
#' Note: I may change this to a ggplot and provide options to control sliding windows
#' 
#' @param dna Character string, or a DNAString object, or a length 1 DNAStringSet object as defined by the Biostrings package.
#' @param rp Text string defining that defines how complexity is defined. [AT PRESENT ONLY PATTERNS COMPOSED OF "1" and "." CAN BE USED]. Default is "1.1", which specifies to return the fraction of bases in which base[i] == base[i+2].
#' @param wsize Number that controls windowsizes; default window sizes of the larger windows (indicated by red and blue lines) are multiplied by this number.
#' @return NULL
#' @export repeatScorePlot
repeatScorePlot <- function(dna,rp="1.1",wsize=1){
	dna   <- DNAStringSet(dna)
	smallwindow <- min( max((width(dna)/100),nchar(rp)),10)
	smallshift  <- 1
	bigwindow   <- min(max((width(dna)/10),nchar(rp)),50)
	bigshift    <- bigwindow
	test1 <- wapply(x=dna, FUN=function(x){repeatScore(x,rp=rp)}, windowsize=smallwindow, windowshift=smallshift)
	test2 <- wapply(x=dna, FUN=function(x){repeatScore(x,rp=rp)}, windowsize=bigwindow*wsize, windowshift=(bigshift*wsize))
	sim1  <- DNAStringSet(paste(sample(c("A","C","G","T"),size=width(dna),replace=T),collapse=""))
	test3 <- wapply(x=sim1, FUN=function(x){repeatScore(x,rp=rp)}, windowsize=bigwindow*wsize, windowshift=(bigshift*wsize))
	plot(seq(from=1,to=width(dna),length.out=length(test1)),test1,type="l",xlab="nucleotide position along sequence",ylab=sprintf("DNA repeat score for pattern '%s'",rp),sub="black and red = empirical; blue = random",ylim=c(0,1))
	lines(seq(from=1,to=width(dna),length.out=length(test2)),test2,col="red")
	lines(seq(from=1,to=width(sim1),length.out=length(test3)),test3,col="blue")
	res <- list(scores.smallwindow=test1,scores.bigwindow=test2,sim.scores=test3)
	res
}










