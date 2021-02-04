#' Read Multiple Fasta Files
#' 
#' Read fasta formatted sequences stored in multiple files in the same directory.
#' 
#' @param directory Path to folder containing multiple fasta formatted sequence files.
#' @return A DNAStringSet, RNAStringSet, or AAStringSet (determined by the sequences).
#' @export
scanFa.multiple <- function(directory){
	scan.fa <- function(filename){
		res <- Rsamtools::scanFa(Rsamtools::FaFile(filename))
		res
	}
	
	files <- list.files(directory,full.names=T)
	for(i in 1:length(files)){
		temp.data <- Rsamtools::scan.fa(files[i])
		if(i==1){
			all.data <- temp.data
		}
		if(i>1){
			all.data <- append(all.data,temp.data)
		}
	}
	all.data
}
