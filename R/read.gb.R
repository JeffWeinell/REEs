#' Read GenBank flatfile
#' 
#' This function imports a GenBank flatfile into R
#' 
#' @param gb.filename Path to the name GenBank flatfile.
#' @param progress Show a progress bar while importing data into R (default is TRUE).
#' @return Object of class gbRecord that contains the information from the GenBank flatfile.
#' @export 
read.gb <- function(gb.filename,progress=T){	
	include.accession.region.gb <- function(gb.filename, outfile){
		original.string  <- readChar(gb.filename, file.info(gb.filename)$size)
		new.string       <- gsub(pattern=" REGION: ",replacement=":",original.string)
		write(new.string,file=outfile)
	}
	if(progress==T){
		progress.val <- T
	} else {
		progress.val <- F
	}
	outfile.tmp    <- paste(gsub("\\..+","",gb.filename),".tmp",sep="")
	include.accession.region.gb(gb.filename, outfile=outfile.tmp)
	gbData         <- biofiles::gbRecord(outfile.tmp,progress=progress.val)   ### read in the GenBank flatfile (which can contain multiple records)
	if(gb.filename!=outfile.tmp){
		file.remove(outfile.tmp)
	}
	gbData
}
