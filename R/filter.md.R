#' @title Test Missing Data
#'
#' Test if a character vector contains. This function is used internally in the the filter.alignment function.
#'
#' @param alignment.site Character vector corresponding to a site in a DNA or AA alignment
#' @param missing.data.threshold Value between 0 and 1 and determines the maximum fraction of missing characters in the input vector for the function to return TRUE.
#' @param type Either "DNA" or "AA"
#' @param ambiguous Logical indicating if ambiguous characters should be treated as missing data (default is FALSE)
#' @return Returns TRUE if no missing data, or if the fraction of missing data is less than the value of missing.data.threshold argument; otherwise returns FALSE.
#' @export filter.md
filter.md <- function(alignment.site,missing.data.threshold,type,ambiguous=F) {
	## This is the number of each type of character in a given column
	#y              <- table(alignment.site)
	## Number of each type of character in a given column
	allele.counts   <- table(tolower(alignment.site))
	
	if(ambiguous==F){
		fraction.missing <- (sum(allele.counts[which(names(allele.counts) %in% c("-"))])/sum(allele.counts))
		md.y <- fraction.missing
	}
	if(type=="DNA" & ambiguous==T){
		fraction.missing.ambiguous <- (1-(sum(allele.counts[(names(allele.counts) %in%  c("a","c","g","t"))]) /sum(allele.counts)))
		md.y <- fraction.missing.ambiguous
	}
	if(type=="AA" & ambiguous==T){
		fraction.missing.ambiguous <- (1-(sum(allele.counts[(names(allele.counts) %in% tolower(REEs::s2v("*ACDEFGHIKLMNPQRSTVWY")))])/sum(allele.counts)))
		md.y <- fraction.missing.ambiguous
	}
	
	if(md.y < missing.data.threshold | md.y==0){
		z = TRUE
	} else {
		z = FALSE
	}
	z
}
