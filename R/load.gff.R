#' Read a General Feature Format (GFF) file
#' 
#' This function reads a General Feature Format (GFF) file. If the GFF is a local file, then this function is simply a wrapper for the function read.gff() in the package ape.
#' If the GFF is stored in an open source online data repository (e.g. NCBI) this function first downloads a temporary local copy of the file.
#' See http://gmod.org/wiki/GFF3 for info about GFF files.
#' 
#' @param input A character string giving the path to the GFF file. If input is a URL then local must be FALSE.
#' @param local Whether or not input is a path to a local file (default is TRUE). If set to FALSE, then the function treats input as URL.
#' @return A data table object containing the GFF feature rows. Incomplete feature rows are dropped.
#' @export load.gff
load.gff <- function(input,local=TRUE){
	if(local){
		gff.obj <- ape::read.gff(input)
	} else {
		### Creates a path (character string) that the GFF will be saved to 
		temp.name <- tempfile()
		### Downloads the GFF to temp.name
		utils::download.file(input,dest=temp.name,quiet=T)
		# Reads the GFF into R
		gff.obj   <- ape::read.gff(temp.name, na.strings = c(".", "?"), GFF3 = TRUE)
		# Removes the temporary file
		invisible(file.remove(temp.name))
	}
	gff.obj <- data.table::as.data.table(gff.obj)
	gff.obj
}

#' @examples
#' Thamnophis.sirtalis_GFF.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"
#' Thamnophis.sirtalis_GFF     <- REEs::load.gff(input=Thamnophis.sirtalis_GFF.url,local=F)
#' 
