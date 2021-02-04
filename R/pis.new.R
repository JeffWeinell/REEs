#' Parsimony Informative Sites
#' 
#' Calculates either the number of parsimony informative sites or the fraction or percent of sites that parsimony inforative.
#' 
#' @param x input alignment of DNA sequences stored as either a DNAStringSet, DNAMultipleAlignment, or DNAbin. This is an improved version of the pis function from the phyloch package.
#' @param abs Should the output of this function be the absolute number of parsimony informative sites? Default is TRUE. If FALSE, the argument as.percent determines the output type.
#' @param use.ambiguities Should ambiguous bases be considered as non-missing data? Default is false (Don't change this).
#' @param as.percent When abs = FALSE, should the output of this function be a percent? Default is FALSE, in which case the output is the fraction of sites parsimony informative.
#' @return If abs = TRUE, this function returns the number of parsimony informative sites in the alignment.
#'  If abs=FALSE and as.percent = TRUE, this function returns the percent of sites that are parsimony informative.
#'  If abs=FALSE and as.percent = FALSE, this function returns the fraction of sites that are parsimony informative.
#' @export 
pis.new <- function (x, abs = TRUE, use.ambiguities = FALSE,as.percent=F){
	if(class(x) %in% c("DNAStringSet","DNAMultipleAlignment")){
		x <- ape::as.DNAbin(Biostrings::DNAMultipleAlignment(x))
	}
	nbchar <- dim(x)[2]
	x   <- as.character(x)
	out <- apply(x, 2, REEs::pars.inf)
	out <- length(out[out])
	
	if (!abs){
		out <- round(out/nbchar * 100, digits = 2) #### Output is the percentage of sites with at least two different and not ambiguous character states
	}
	if(abs ==F & as.percent==F){
		out <- out/100
	}
	out
}

#' Parsimony Informative Sites
#' 
#' Calculates either the number of parsimony informative sites or the fraction or percent of sites that parsimony inforative.
#' This is the same function as used in the ips package.
#' 
#' @param x input alignment of DNA sequences stored as either a DNAStringSet, DNAMultipleAlignment, or DNAbin. This is an improved version of the pis function from the phyloch package.
#' @param what Either of "absolute", "fraction", or "index", which will return the absolute number, the relative number or the indeces of the potentially-informative sites.
#' @param use.ambiguities Should ambiguous bases be considered as non-missing data? Default is false (Don't change this).
#' @return Numeric (depending on what, the number, fraction, or indices of potentially-informative nucleotide sites).
#' @export pis
pis <- function (x, what = "fraction", use.ambiguities = FALSE){
	if(!inherits(x, "DNAbin")){
		stop("'x' is not of class 'DNAbin'")
	}
	what <- match.arg(what, c("absolute", "fraction", "index"))
	if(use.ambiguities){
		warning("'use.ambiguities' is currently ignored ", "and IUPAC ambiguity symbols are treated as missing data")
		use.ambiguities <- FALSE
	}
	x   <- as.character(x)
	out <- apply(x, 2, REEs::pars.inf)
	if (what %in% c("absolute", "fraction")){
		out <- length(out[out])
		if (what == "fraction") {
			out <- round(out/ncol(x) * 100, digits = 2)
		}
	}
	else {
		out <- which(out)
	}
	out
}


#' Is Site Parsimony Informative
#' 
#' Given a vector of homologous DNA characters, returns TRUE if the site is parsimony informative and FALSE otherwise.
#' This is the pars.inf function that was used internally in the pis function of the ips package.
#' Added "?" to the list of possible characters treated as missing data.
#' 
#' @param x Character vector of DNA at one site of an alignment.
#' @return TRUE if the site is parsimony informative. Otherwise, FALSE.
#' @export pars.inf
pars.inf <- function(x) {
		x <- table(x)
		# This is the number of each type of character in a given column
		x <- x[x > 1]
		# A list of ambiguous or missing characters to check against.
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")
		if (length(x[!names(x) %in% n]) > 1){
			x <- TRUE
		} else {
			x <- FALSE
		}
	}

