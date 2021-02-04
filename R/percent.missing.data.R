#' Percent Missing Data
#' 
#' Calculates the percent of bases (not sites) in a DNA alignment that are missing (i.e., gaps or ambiguous).
#' 
#' @param alignment DNA alignment stored as an DNAStringSet
#' @return Percent of data missing (gaps or ambiguous) in a DNA sequence alignment
#' @export 
percent.missing.data <- function(alignment){
	store.class <- class(alignment)
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- as.DNAbin(DNAMultipleAlignment(alignment))
		data.type <- "DNA"
	}
	if(store.class %in% c("AAMultipleAlignment","AAStringSet")){
		alignment <- as.AAbin(AAMultipleAlignment(alignment))
		data.type <- "AA"
	}
	x      <- as.character(alignment)
	y      <- table(x)
	n      <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")
	y.n    <- y[(names(y) %in% n)]
	result <- round((sum(y.n)/sum(y)*100),digits=2)
	result
}
