#' Trim DNA alignment
#' 
#' Trims the first and last columns of a DNA alignment until those columns contain at least the defined threshold number of individuals with non-missing data.
#' 
#' @param  alignment DNA alignment held in an object of class DNAStringSet, DNAMultipleAlignment, or DNAbin. These object classes are defined in the Biostrings and ape packages.
#' @param  threshold Number of individuals with non-missing data that must be present to keep the alignment ends. Default is 1.
#' @return Trimmed alignment object having the same class as the input alignment.
#' @export 
trim.alignment <- function (alignment, threshold=1){
	
	store.class <- class(alignment) ### Object class of the input alignment. This is used so that the output alignment can have the same class as the input alignment.
	
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- ape::as.DNAbin(Biostrings::DNAMultipleAlignment(alignment))
	}
	
	# Tests if an alignment column has greater than or equal to x (=threshold) number of non-missing data (A,C,G,and Ts)
	meets.threshold <- function(x) {
		x <- table(x)                                                     # This is the number of each type of character in a given column
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")  # A list of possible characters to check against
		if (sum(x[!names(x) %in% n]) > threshold){                        #| These lines are specifying whether the 
			x <- TRUE                                                     #| the ith column has at least theshold characters that are not ambiguous or gaps
		} else {                                                          #|
			x <- FALSE                                                    #| the ith column has fewer than the theshold number of non-ambiguous/non-gap characters
		}                                                                 #|
	}
	
	# Commented out because this function is loaded in as stand-alone function, so no need to include it as an internal function.
	#string.list <- function(y){
	#	paste.collapse <- function(z){
	#		z <- paste(z,collapse="")
	#		z
	#	}
	#	y <- as.character(y)
	#	y <- apply(X=y,MARGIN=1,FUN=paste.collapse)
	#	y
	#}
	
	x         <- as.character(alignment)
	logicals  <- unlist(apply(x, 2, meets.threshold))
	new.start <- min(which(logicals))
	new.end   <- max(which(logicals))
	out       <- alignment[,new.start:new.end]

	keep.rows <- apply(X=out,MARGIN=1,FUN=charTest.dna,method="any") ### A list of logicals indicating whether each individual has at least one non-missing character
	out       <- out[keep.rows,]                                     ### Only keeps individual that have some data present
	if(store.class=="DNAStringSet"){
		out   <- Biostrings::DNAStringSet(string.list(out))
	}
	if(store.class=="DNAMultipleAlignment"){
		out   <- Biostrings::DNAMultipleAlignment(string.list(out))
	}
	out
}
