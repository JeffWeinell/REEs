#' Get SNP Alignment
#' 
#' Filters an input alignment and returns an alignment containing only the sites that are SNPs.
#' This function is the same as running filter.alignment(alignment,min.allele.freqs.dna=c(2,2,0,0))
#' 
#' @param alignment Input alignment
#' @return An alignment containing only the SNPs of the input alignment.
#' @export 
snp.alignment <- function (alignment){
	
	######
	## Defining internal functions needed
	## for the whole function to work
	#### pars.inf.dna function ### tests if an alignment site is parsimony informative
	#pars.inf.dna <- function(alignment.site) {	   # alignment.site is one column of a DNA alignment
	#	y <- table(alignment.site)                 # This is the number of each type of character in a given column
	#	n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")  # A list of possible ambiguous or missing characters to ignore
	#	y <- y[!(names(y) %in% n)]
	#	if(any(y>1) & length(y)>1){
	#		y <- TRUE
	#	} else {
	#		y <- FALSE
	#	}
	#	y
	#}	
	##### pars.inf.aa function
	#pars.inf.aa <- function(alignment.site) {	   #|alignment.site is a column of an AA alignment
	#	y   <- table(alignment.site)               #|this is the number of each type of character in a given column
	#	n   <- c("-","X")                          #|a list of possible ambiguous or missing characters to ignore
	#	AAs <- c(unique(GENETIC_CODE))
	#	y   <- y[!(names(y) %in% n)]
	#	if(any(y>1) & length(y)>1){                #|tests if the column has at least two different characters,
	#		y <- TRUE                              #|and if at least on of the characters if shared by more than one individual
	#	} else {
	#		y <- FALSE
	#	}
	#	y
	#}
	#### string.list function
	#string.list <- function(y){
	#	paste.collapse <- function(z){
	#		z <- paste(z,collapse="")
	#		z
	#	}
	#	if(class(y)!="matrix"){
	#		y <- as.character(y)
	#	}
	#	y <- apply(X=y,MARGIN=1,FUN=paste.collapse)
	#	y
	#}
	#### charTest.dna function ### tests if an individual has at least one of each nucleotide
	#charTest.dna <- function(ls){
	#	ls <- tolower(ls)
	#	any(c("a","c","g","t") %in% ls)  ### changed "all" to "any"
	#}
	#### charTest.aa function
	#charTest.aa <- function(ls){
	#	#ls <- tolower(ls)
	#	#all(c("a","c","g","t") %in% ls)
	#	any(GENETIC_CODE %in% ls)
	#}
	
	#### Start of actual function  ####
	store.class <- class(alignment)
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- ape::as.DNAbin(Biostrings::DNAMultipleAlignment(alignment))
		data.type <- "DNA"
	}
	if(store.class %in% c("AAMultipleAlignment","AAStringSet")){
		alignment <- ape::as.AAbin(Biostrings::AAMultipleAlignment(alignment))
		data.type <- "AA"
	}
	x <- as.character(alignment)
	if(data.type == "DNA"){
		keep.cols  <- unlist(apply(X=x, MARGIN=2, FUN=pars.inf.dna))            ### Tests if a site is parsimony informative
		out        <- x[,keep.cols]
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.dna,method="any") ### A list of logicals indicating whether each individual has at least one non-missing character

	}
	if(data.type == "AA"){
		keep.cols  <- unlist(apply(X=x, MARGIN=2, FUN=pars.inf.aa))          ### Tests if a site is parsimony informative
		out        <- x[,keep.cols]
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.aa)     ### A list of logicals indicating whether each individual has at least one non-missing character
	}
	out        <- out[keep.rows,]                               ### Only keeps individuals that have some data present
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

#' SNPs DNAStringSet Alignment
#' 
#' Removes non-SNP columns of a DNA alignment of class DNAStringSet. Does the same thing as snp.alignment function, but this function is more efficient although only works on DNAStringSet objects.
#' 
#'  @param alignment DNA sequence alignment of class DNAStringSet
#'  @param remove.NoDataSeqs Filter the SNPs alignment to only include individuals with at least one A, C, G, or T.
#'  @return DNAStringSet alignment of SNPs
#'  @export 
snpsDNAStringSet <- function(alignment,remove.NoDataSeqs = T){
	keep.cols <- rep(F,width(alignment[1]))
	for(i in 1:width(alignment[1])){
		chars.temp <- table(Biostrings::subseq(alignment,start=i,width=1))
		nucs.temp  <- chars.temp[names(chars.temp) %in% c("A","C","G","T") & chars.temp > 1]
		if(length(nucs.temp)>1){
			keep.cols[i] <- T
		}
	}
	alignment.mat       <- do.call(rbind,strsplit(as.character(alignment),split=""))
	alignment.mat.snps  <- alignment.mat[,keep.cols]
	alignment.list.snps <- apply(X=alignment.mat.snps,MARGIN=1,FUN=function(y){paste(y,collapse="")})
	alignment.set.snps.temp  <- Biostrings::DNAStringSet(alignment.list.snps)
	if(remove.NoDataSeqs == T){
		keep.rows          <- which(countChars(alignment.set.snps.temp,"A|C|G|T")>0)
		alignment.set.snps <- alignment.set.snps.temp[keep.rows]
	} else {
		alignment.set.snps <- alignment.set.snps.temp
	}
	alignment.set.snps
}

