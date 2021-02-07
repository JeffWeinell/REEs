#' Read a General Feature Format (GFF) file
#' 
#' This function reads a General Feature Format (GFF) file. If the GFF is a local file, then this function is simply a wrapper for the function read.gff() in the package ape.
#' If the GFF is stored in an open source online data repository (e.g. NCBI) this function first downloads a temporary local copy of the file.
#' See http://gmod.org/wiki/GFF3 for info about GFF files.
#' 
#' @param input.gff One of the following: An object of class data table, data.frame, or character matrix, with columns identical to those used in Generic Feature Format (GFF) tables; or, a character string with local file path to either a GFF file or table with columns identical to those used in GFF files; or, a character string with URL path to a GFF file.
#' @return A data table object containing the GFF feature rows. Incomplete feature rows are dropped.
#' @export load.gff
load.gff <- function(input.gff){
	if(is(input.gff,"character")){
		if(file.exists(input.gff)){
			### Check first line of file to determine if it is a GFF-like table or an actual GFF file (with a GFF header)
			Line1 <- readLines(input.gff,n=1)
			if(any(grep("gff-version",Line1))){
				### Read GFF using ape function read.gff
				gff.obj <- ape::read.gff(input.gff, na.strings = c(".", "?"), GFF3 = TRUE)
			} else {
				### Read the GFF as a data.table object and then coerce it to a data frame.
				gff.obj <- data.table::fread(input.gff)
			}
		} else {
			### Creates a path (character string) that the GFF will be saved to 
			temp.name <- tempfile()
			### Downloads the GFF to temp.name
			utils::download.file(input.gff,dest=temp.name,quiet=T)
			### Now check if the downloaded file is a GFF-like table or an actual GFF file (with a GFF header)
			Line1 <- readLines(temp.name,n=1)
			if(any(grep("gff-version",Line1))){
				# Reads the GFF into R as a data.frame object
				gff.obj   <- ape::read.gff(temp.name, na.strings = c(".", "?"), GFF3 = TRUE)
				# Removes the temporary file
				invisible(file.remove(temp.name))
			} else {
				gff.obj <- data.table::fread(temp.name)
			}
		}
	} else {
		if(is(input.gff,"data.frame") | is(input.gff,"matrix")){
			gff.obj   <- input.gff
		} else {
			stop("input.gff has unrecognized format")
		}
	}
	### Coerce to data.table object
	gff.obj <- data.table::as.data.table(gff.obj)
	gff.obj
}
#' @examples
#' Thamnophis.sirtalis_GFF.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"
#' Thamnophis.sirtalis_GFF     <- REEs::load.gff(input.gff=Thamnophis.sirtalis_GFF.url,local=F)
#' 
##### Previous version of function is shown below.
# @param local Whether or not input.gff is a path to a local file (default is TRUE). If set to FALSE, then the function treats input.gff as URL.
#########
#	if(local){
#		gff.obj <- ape::read.gff(input.gff)
#	} else {
#		### Creates a path (character string) that the GFF will be saved to 
#		temp.name <- tempfile()
#		### Downloads the GFF to temp.name
#		utils::download.file(input.gff,dest=temp.name,quiet=T)
#		# Reads the GFF into R
#		gff.obj   <- ape::read.gff(temp.name, na.strings = c(".", "?"), GFF3 = TRUE)
#		# Removes the temporary file
#		invisible(file.remove(temp.name))
#	}
#	gff.obj <- data.table::as.data.table(gff.obj)
#	gff.obj
#}

