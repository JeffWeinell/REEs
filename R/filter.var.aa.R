#' Filter by AA Variation
#' 
#' This function is used internally in the function filter.alignment to test if an alignment site meets threshold levels of amino acid variation.
#' Default values return TRUE for any input character vector.
#' 
#' @param alignment.site Input character vector representing AA data at one site of an alignment.
#' @param min.freqs Numerical vector. Values determine the minimum frequency of the most frequent, second-most frequent, ..., and n-most frequent allele type (amino acid) required to return a value of TRUE.
#'  Additionally: (1) the length of min.freqs determines the maximum number of different AAs allowed or else function returns FALSE.
#'                (2) the number of non-zero elements in min.freqs determines the minimum number of different AAs required or else the function returns FALSE.
#'  Default value is a vector of zeros with length 21, meaning that up to 21 allele types (the 20 amino acids and stop codon) are allowed, and that none of these AAs need to be present for the function to return a value of TRUE.
#' @return Logical indicating if the site meets all of the input requirements for amino acid variation.
#' @export filter.var.aa
filter.var.aa <- function(alignment.site,min.freqs=rep(0,21)){
	allele.counts  <- table(tolower(alignment.site))        ## Number of each type of character in a given column
	AA.counts      <- allele.counts[(names(allele.counts) %in% tolower(REEs::s2v("*ACDEFGHIKLMNPQRSTVWY")))]
	AA.counts.sorted <- sort(AA.counts,decreasing=T)
	if(length(AA.counts) <= length(min.freqs) & length(AA.counts) >= length(which(min.freqs > 0))){
		if(length(which(min.freqs > 0))==0){
			result <- TRUE
		} else {
			if(all((AA.counts.sorted[min.freqs>0]-min.freqs[min.freqs>0])>=0)){
				result <- TRUE
			} else {
				result <- FALSE
			}
		}
	} else {
		result <- FALSE
	}
	result
}


#' Test if AA Site Parsimony Informative
#' 
#' This function is a special case of the function filter.var.aa and is used internally in the function snp.alignment to test if an amino acid alignment site is parsimony informative.
#' Default values return TRUE for any input character vector.
#' 
#' @param alignment.site Input character vector representing AA data at one site of an alignment.
#' @return Logical indicating if the site is parsimony informative.
#' @export pars.inf.aa
pars.inf.aa <- function(alignment.site){
	input  <- alignment.site
	result <- filter.var.aa(input,min.freqs=c(c(2,2),rep(0,19)))
	result
}
