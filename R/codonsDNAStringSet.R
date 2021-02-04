#' Get Sites at Codon Position
#' 
#' Returns a DNA alignment of only the first, second, or third codon position sites (or any two of these three).
#' 
#' @param x Input alignment DNAStringSet object.
#' @param which.codon.positions An integer in {1,2,3}, or an integer vector containing the values 1, 2, and/or 3.
#' @param predict.reading.frame Whether or not to try to guess the reading frame, based on the relative number of snps in each frame. Default is FALSE.
#'  The third frame is assumed to be the frame with the most snps. This method works best for longer, faster-evolving loci, and when the alignment has many individuals, and/or moderate to high genetic distances among individuals.
#' @param frameAtFirstchar The frame (i.e., codon position) of the first character of the sequence (if known). If unknown, leave as NULL, in which case the first base treated as if it is the first codon position.
#' @return DNAStringSet object holding only the sites at codon positions indicated by the argument which.codon.positions
#' @export 
codonsDNAStringSet <- function(x,which.codon.positions,predict.reading.frame=F,frameAtFirstchar=NULL){
	x.mat               <- do.call(rbind,strsplit(as.character(x),split=""))
	if(predict.reading.frame==T){
		x.frame1.mat    <- x.mat[,c(T,F,F)]
		x.frame2.mat    <- x.mat[,c(F,T,F)]
		x.frame3.mat    <- x.mat[,c(F,F,T)]
		x.frame1.set    <- Biostrings::DNAStringSet(apply(X=x.frame1.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
		x.frame2.set    <- Biostrings::DNAStringSet(apply(X=x.frame2.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
		x.frame3.set    <- Biostrings::DNAStringSet(apply(X=x.frame3.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
		snps.per.frame  <- c(width(snpsDNAStringSet(x.frame1.set)[1]),width(snpsDNAStringSet(x.frame2.set)[1]),width(snpsDNAStringSet(x.frame3.set)[1]))
		frame3          <- which(snps.per.frame==max(snps.per.frame))
		frame1          <- c(frame3-2,frame3+1)[c(frame3-2,frame3+1) %in% c(1:3)]
		frame2          <- c(frame3-1,frame3+2)[c(frame3-1,frame3+2) %in% c(1:3)]
		frames          <- c(frame1,frame2,frame3)
		codon.positions <- frames %in% which.codon.positions
	}
	if(is.null(frameAtFirstchar)==F){
		char1.frame     <- frameAtFirstchar
		char2.frame     <- c(char1.frame-2,char1.frame+1)[c(char1.frame-2,char1.frame+1) %in% c(1:3)]
		char3.frame     <- c(char2.frame-2,char2.frame+1)[c(char2.frame-2,char2.frame+1) %in% c(1:3)]
		frames          <- c(char1.frame,char2.frame,char3.frame)
		codon.positions <- frames %in% which.codon.positions
	}
	if(predict.reading.frame==F & is.null(frameAtFirstchar)==T) {
		codon.positions <- c(1:3) %in% which.codon.positions
	}
	positions.mat <- x.mat[,codon.positions]
	positions.set <- Biostrings::DNAStringSet(apply(X=positions.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
	positions.set
}
