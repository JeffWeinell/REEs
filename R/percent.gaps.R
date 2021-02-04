#' Percentage of Alignment Gaps
#' 
#' Calculates the percentage of bases (not sites) in a DNA alignment that are gaps.
#' 
#' @param alignment DNA alignment stored as an DNAStringSet
#' @return Percentage of bases that are gaps
#' @export 
percent.gaps <- function(alignment){
	alignment <- Biostrings::as.DNAStringSet(alignment)
	area.alignment <- ((unique(width(alignment)))*length(alignment))
	### Pass the alignment to C++ compiled DECIPHER function to remove gaps quickly.
	noGaps.seqs    <- DECIPHER::RemoveGaps(alignment)
	result         <- (area.alignment-sum(width(noGaps.seqs)))/area.alignment
	result
}