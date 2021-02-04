#' IUPAC nucleotide symbols
#' 
#' This function returns the IUPAC symbol for a nucleotide sequence, for instance c("c", "c", "g") is coded by "s",
#' which is useful for describing ambiguous or consensus nucleotide sequences.
#' This function is identical to the function bma in package seqinr.
#' 
#' @param nucl A nucleotide sequence as a vector of single characters (e.g., nucl=c("c", "c", "g")).
#' @param warn.non.IUPAC If TRUE (default) warns when no IUPAC symbol is possible.
#' @param type Type of nucleotides: either "DNA" or "RNA"
#' @return A single IUPAC symbol in lower case, or NA when this is not possible.
#' @export bma
bma <- function (nucl, warn.non.IUPAC = TRUE, type = c("DNA", "RNA")){
	### Check if nucl is a vector of single characters and stop with warning if not.
	if (nchar(nucl[1]) != 1){
		stop("vector of single chars expected")
	}
	### Takes the first entry in the value of type, which means that DNA is the default of the type argument.
	type  <- match.arg(type)
	### Store nucleotide characters as lower case
	nucl  <- tolower(nucl)
	nucl  <- unlist(sapply(X=nucl, FUN=REEs::amb, checkBase = FALSE))
	iupac <- sapply(REEs::amb(), REEs::amb)
	if (type == "DNA") {
		iupac$u <- NULL
	} else {
		iupac$t <- NULL
	}
	idx <- unlist(lapply(X=iupac, FUN=setequal, y= nucl))
	if (all(idx == FALSE)) {
		if (warn.non.IUPAC) {
			warning(paste("Undefined IUPAC symbol with:", paste(nucl, collapse = " ")))
		}
		return(NA)
	}
	return(names(iupac)[idx])
}

#' Expansion of IUPAC nucleotide symbols
#' 
#' This function returns the list of nucleotide matching a given IUPAC nucleotide symbol, for instance c("c", "g") for "s".
#' The function is a slightly modified version (internal syntax only) of the function amb in the package seqinr.
#' 
#' 
#' @param base An IUPAC symbol for a nucleotide as a single character
#' @param forceToLower If TRUE (the default) the base is forced to lower case.
#' @param checkBase If TRUE (the default) the character is checked to belong to the allowed IUPAC symbol list.
#' @param IUPAC A character vector containing the allowed IUPAC symbols.
#' @param u2t If TRUE (the default), "u" for uracil in RNA are changed into "t" for thymine in DNA.
#' @return When base is missing, the list of IUPAC symbols is returned, otherwise a vector with expanded symbols.
#' @export amb
amb <- function (base, forceToLower = TRUE, checkBase = TRUE, IUPAC = REEs::s2v("acgturymkswbdhvn"), u2t = TRUE){
	# IUPAC <- REEs::s2v("acgturymkswbdhvn")
	#unlist(strsplit("acgturymkswbdhvn",split=""))
	### Chck if base has a value. If not, set base equal to IUPAC.
	if (missing(base)) 
		return(IUPAC)
	if (!is.character(base)){
		stop("Character expected")
	}
	if (nchar(base) != 1){
		stop("Single character expected")
	}
	if (forceToLower){
		base <- tolower(base)
	}
	if (checkBase){
		if (!(base %in% IUPAC)){
			stop("IUPAC base expected")
		}
	}
	bases            <- c("r","y","m","k","s","w","b","d","h","v","n","u")
	if(u2t){
		names(bases) <- c("ag","ct","ac","gt","cg","at","cgt","agt","act","acg","acgt","t")
	} else {
		names(bases) <- c("ag","ct","ac","gt","cg","at","cgt","agt","act","acg","acgt","u")
	}
	return(REEs::s2v(names(bases[which(bases %in% base)])))
}

#' Character String to Character Vector
#' 
#' Split all characters in a character string into a vector of single characters.
#' The is lightly modified version (internal syntax only) of the function s2c function in the package seqinr.
#' @param string One character string.
#' @return A character vector with the ith entry equal to the ith character in the input string. If length(string) is not zero or if string contains non-characters, this function returns a warning and NA.
s2v <- function(string){
	if((!is(string,"character")) | (length(string)!=1)){
		warning("string must but a character vector of length 1, returning NA")
		return(NA)
	} else {
		return(unlist(strsplit(string,split="")))
	}
}

