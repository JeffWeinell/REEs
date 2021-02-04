#' Filter by DNA Variation
#' 
#' This function is used internally in the function filter.alignment to test if a DNA alignment site is meets threshold levels of variation.
#' Default values return TRUE for any input character vector.
#' 
#' @param alignment.site Input character vector representing DNA data at one site of an alignment.
#' @param min.freqs Numerical vector of length n, indicating the minimum frequencies of the most frequent, second-most frequent, third-most frequent, and fourth-most frequent allele required to return a value of TRUE.
#'  Additionally, the length of min.freqs sets the maximum number of allele types required, and the number of non-zero elements in min.freqs sets the minimum number of allele types required.
#'  Default is c(0,0,0,0), meaning that up to four allele types are allowed (the four DNA nucleotides), and that none of these need to be present for the function to return a value of TRUE.
#' @return Logical indicating if the site meets all of the input requirements for nucleotide variation.
#' @export filter.var.dna
filter.var.dna <- function(alignment.site,min.freqs=c(0,0,0,0)){
	allele.counts   <- table(tolower(alignment.site))        ## Number of each type of character in a given column
	ACGT.counts     <- allele.counts[(names(allele.counts) %in% c("a","c","g","t"))]
	ACGT.counts.sorted <- sort(ACGT.counts,decreasing=T)
	if(length(ACGT.counts) <= length(min.freqs) & length(ACGT.counts) >= length(which(min.freqs > 0))){
		if(length(which(min.freqs > 0))==0){
			result <- TRUE
		} else {
			if(all((ACGT.counts.sorted[min.freqs>0]-min.freqs[min.freqs>0])>=0)){
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



#' Test if DNA Site Parsimony Informative
#' 
#' This function is a special case of the function filter.var.dna and is used internally in the function snp.alignment to test if a DNA alignment site is parsimony informative.
#' Default values return TRUE for any input character vector.
#' 
#' @param alignment.site Input character vector representing DNA data at one site of an alignment.
#' @return Logical indicating if the site is parsimony informative.
#' @export pars.inf.dna
pars.inf.dna <- function(alignment.site){
	input  <- alignment.site
	result <- filter.var.dna(input,min.freqs=c(2,2,0,0))
	result
}


