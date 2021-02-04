#' Filter Alignment
#' 
#' Filters individuals or sites from an DNA or AA alignment based on user-defined thresholds of percent missing or ambiguous data per site, variation per site, frequency of common alleles per site, frequency of rare alleles per site, or specific sites (e.g., every second, third, etc.)
#' Running this function with the default values for arguments returns the input alignment unchanged.
#' 
#' @param alignment Input alignment of class DNAStringSet, DNAMultipleAlignment, or DNAbin
#' @param mdt Missing data threshold: maximum fraction of missing data allowed at a site (default = 1)
#' @param treat.ambiguous.as.missing Logical indicating if ambiguous characters should be treated as missing data (default is FALSE).
#' @param min.allele.freqs.dna Numerical vector of length n, indicating the minimum frequencies of the most frequent, second-most frequent, third-most frequent, and fourth-most frequent allele required to return a value of TRUE.
#'  Additionally, the length of the input vector sets the maximum number of allele types allowed. The number of non-zero values sets the minimum number of allele types required.
#'  Default is c(0,0,0,0), meaning that up to four different DNA nucleotides are allowed, and that none of these need to be present for the function to return a value of TRUE for the site.
#' @param min.allele.freqs.aa Numerical vector. Values determine the minimum frequency of the most frequent, second-most frequent, ..., and n-most frequent allele type (amino acid) required to return a value of TRUE.
#'  Additionally: (1) the length of min.freqs determines the maximum number of different AAs allowed or else function returns FALSE.
#'                (2) the number of non-zero elements in min.freqs determines the minimum number of different AAs required or else the function returns FALSE.
#'  Default value is a vector of zeros with length 21, meaning that up to 21 allele types (the 20 amino acids and stop codon) are allowed, and that none of these AAs need to be present for the function to return a value of TRUE.
#' @param keep.specific Providing a numeric range or set allows you to extract specific sites (default NA).
#' @param drop.nodata.samples Drops an individual if it doesnt have any data for the sites that meet the other criteria (default FALSE).
#' @return A filtered DNA or AA alignment
#' @export filter.alignment
filter.alignment <- function(alignment,mdt=1,treat.ambiguous.as.missing=F,min.allele.freqs.dna=c(0,0,0,0),min.allele.freqs.aa=rep(0,21),keep.specific=NA,drop.nodata.samples = F){

	store.class <- class(alignment)
	
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- ape::as.DNAbin(Biostrings::DNAMultipleAlignment(alignment))
	}
	
	if(store.class %in% c("AAMultipleAlignment","AAStringSet")){
		alignment <- ape::as.AAbin(Biostrings::AAMultipleAlignment(alignment))
	}
	
	if(class(alignment) == "DNAbin"){
		data.type <- "DNA"
	}
	
	if(class(alignment) == "AAbin"){
		data.type <- "AA"
	}
	
	x <- as.character(alignment)
	
	if(is.na(all(keep.specific))){
		keep.cols3 <- rep(TRUE,ncol(x))
	} else {
		keep.cols3 <- c(1:ncol(x)) %in% keep.specific
	}

	if(data.type == "DNA"){
		#keep.cols1 <- unlist(apply(X=x, MARGIN=2, FUN=pars.inf.dna))
		keep.cols1 <- unlist(apply(X=x, MARGIN=2, FUN=filter.var.dna,min.freqs=min.allele.freqs.dna))
		keep.cols2 <- unlist(apply(X=x, MARGIN=2, FUN=filter.md,type="DNA",missing.data.threshold=mdt,ambiguous=treat.ambiguous.as.missing))
		keep.cols  <- keep.cols1 & keep.cols2 & keep.cols3                                      ### TRUE if keep.cols1, keep.cols2, and keep.cols3 are all TRUE for a particular site
		out        <- as.matrix(x[,keep.cols])
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.dna,method="any") ### A list of logicals indicating whether each individual has at least one non-missing and non-ambiguous character
	}
	
	if(data.type == "AA"){
		#keep.cols1 <- unlist(apply(X=x, MARGIN=2, FUN=pars.inf.aa))
		keep.cols1 <- unlist(apply(X=x, MARGIN=2, FUN=filter.var.aa,min.freqs=min.allele.freqs.aa))
		keep.cols2 <- unlist(apply(X=x, MARGIN=2, FUN=filter.md,type="AA",missing.data.threshold=mdt,ambiguous=treat.ambiguous.as.missing))
		keep.cols  <- keep.cols1 & keep.cols2 & keep.cols3     ### TRUE if keep.cols1, keep.cols2, and keep.cols3 are all TRUE for a particular site
		out        <- as.matrix(x[,keep.cols])
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.aa)    ### A list of logicals indicating whether each individual has at least one non-missing and non-ambiguous character
	}

	if(drop.nodata.samples == T){
		out        <- out[keep.rows,]                          ### Only keep individuals that have some data present (i.e., at least one A, C, G, or T).
	}
	
	if(store.class=="DNAbin"){
		out <- ape::as.DNAbin(out)
	}
	if(store.class=="DNAStringSet"){
		out <- Biostrings::DNAStringSet(string.list(out))
	}
	if(store.class=="DNAMultipleAlignment"){
		out <- Biostrings::DNAMultipleAlignment(string.list(out))
	}
	if(store.class=="AAbin"){
		out <- ape::as.AAbin(out)
	}
	if(store.class=="AAStringSet"){
		out <- Biostrings::AAStringSet(string.list(out))
	}
	if(store.class=="AAMultipleAlignment"){
		out <- Biostrings::AAMultipleAlignment(string.list(out))
	}
	out
}
#' @examples
#' test.alignment <- alignment.DNAStringSet1
#' Returns an alignment identical to the input alignment
#' alignment.noChange <- filter.alignment(test.alignment)
#' 
#' Returns an alignment containing only sites with at least one non-missing or ambiguous character
#' alignment.sitesWithData <- filter.alignment(test.alignment,min.allele.freqs.dna=c(1,0,0,0))
#' 
#' Returns an alignment containing only sites with at least two different non-missing or ambiguous characters (i.e., drops sites that only contain missing or ambiguous characters plus invarient sites)
#' alignment.sitesVariable <- filter.alignment(test.alignment,min.allele.freqs.dna=c(1,1,0,0))
#'
#' Returns an alignment containing only the parsimony informative sites.
#' alignment.pis <- filter.alignment(test.alignment,min.allele.freqs.dna=c(2,2,0,0))
#'
#' Removes individuals if they don't have any data.
#' alignment.NoDataSamplesDropped <- filter.alignment(test.alignment,drop.nodata.samples = T)
#'
#' Removes sites with with more than two DNA nucleotides represented.
#' alignment.RemovedTooVariable <- filter.alignment(test.alignment,min.allele.freqs.dna=c(0,0))
#' 
#' Alignment containing only invarient sites.
#' alignment.InvarientSites <- filter.alignment(test.alignment,min.allele.freqs.dna=c(1))
#' 
#' Filter sites for which the fraction of missing data is above 0.5 (by default, ambiguous characters are not treated as missing data)
#' alignment.mdt.05 <- filter.alignment(test.alignment,mdt=0.5)
#' 
#' Filter sites if there is any missing data (ambiguous characters still allowed)
#' alignment.mdt.0 <- filter.alignment(test.alignment,mdt=0)
#' 
#' Filter sites for which the fraction of missing data (including ambiguous data) is above 0.5
#' alignment.mdt.05 <- filter.alignment(test.alignment,mdt=0.5,treat.ambiguous.as.missing=T)
#' 
#' Filter sites if there is any missing or ambiguous data
#' alignment.mdt.0 <- filter.alignment(test.alignment,mdt=0,treat.ambiguous.as.missing=T)
