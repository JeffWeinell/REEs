#' Add Region Info to GenBank Flatfile
#' 
#' This function is implemented within the read.gb function, and may not have other stand-alone applications.
#' The function creates a temporary file (.tmp) containing a modified version of the input GenBank flatfile, in which the region information (i.e., start and end info)
#' is included as part of the accession name. This modification allows the region info to be acccessed when using the AA.from.gb function.
#' 
#' @param gb.filename 
#' @return A temporary file that is a modified version of the input file that is read into R when running the read.gb function.
#' @export 
include.accession.region.gb <- function(gb.filename){
	original.string         <- readChar(gb.filename, file.info(gb.filename)$size)
	new.string <- gsub(pattern=" REGION: ",replacement=":",original.string)
	outfile.tmp <- gsub(".gb",".tmp",gb.filename)
	write(new.string,file=outfile.tmp)
}
