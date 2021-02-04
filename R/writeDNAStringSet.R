#' Write DNAStringSet
#' 
#' Write a set of DNA sequences to file in sequential fasta or phylip format.
#' 
#' @param x Input DNAStringSet object (see Biostrings package for details on this class)
#' @param filepath Where to save sequences
#' @param append Should the sequences be appended to a file that already exists? (default is FALSE) 
#' @param format Format to save sequences. Available formats are "fasta" (default) or "phylip".
#' @param charsPerLine Maximum number of bp to save per line (default is 100,000 and intended for writing sequences sequentiall rather than interleaved). Only invoked when format is fasta.
#' @return Writes the input DNAStringSet to filepath. Returns the input DNAStringSet x.
#' @export 
writeDNAStringSet <- function(x,filepath,append=F,format="fasta",charsPerLine=100000){
	data         <- ape::as.DNAbin(x)
	param.append <- append
	if(format=="phylip"){
		param.format = "sequential"
	}
	if(format=="fasta"){
		param.format = "fasta"
	}
	ape::write.dna(x=data,file=filepath,format=param.format,append=param.append,nbcol=1,colsep="",colw=charsPerLine)
	x
}
