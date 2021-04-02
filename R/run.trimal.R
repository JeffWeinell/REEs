#' Run trimal
#' 
#' Wrapper for running trimal.
#' "trimAl is a tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment." - (http://trimal.cgenomics.org/).
#' This function takes an input DNA alignment of class DNAStringSet and saves a temporary copy in fasta format. Then, the function runs trimAL on the fasta formatted alignment and saves a temporary trimAL alignment in the current directory.
#' Next, reads the trimAL alignment into R and deletes it from the current directory.
#' Value returned by this function is the trimAL alignment stored as an object of class DNAStringSet.
#' Original version Code written by Carl Hutter?
#' 
#' @param input.align Input DNA alignment (class DNAStringSet)
#' @param method Method to use (default is "auto"). Currently "auto" is the only method implemented in this function, and runs trimal with the "-automated1" flag to "Use a heuristic selection of the automatic method based on similarity statistics." - (http://trimal.cgenomics.org/use_of_the_command_line_trimal_v1.2)
#' @param locus.name Name of the locus held in the input alignment (default is locus.names[i]). If this is kept at the default value then the function should be run within a loop and there needs to be a vector called locus names.
#' @param trimal.exe.dir Directory where the trimal executable is located.
#' @param trimal.exe.name Name of the trimal executable (default is "trimal") 
#' @return trimAL alignment (class DNAStringSet)
#' @export
run.trimal <- function(input.align, method = "auto",locus.name=locus.names[i],trimal.exe.dir,trimal.exe.name="trimal"){
	# input.align  <- rem.align
	# Finds probes that match to two or more contigs
	save.rownames <- names(input.align)                                       ### names of taxa in input alignment
	write.align   <- as.list(as.character(input.align))                       ### a list of character strings (each character string is the alignment DNA sequence of an individual)
	## Writes input alignment as a fasta file in current directory
	locus.name.fa <- paste0(gsub(pattern = "\\..*", "", locus.name), ".fa")   ### name of temportary fasta alignment that is the same as input.align
	seqinr::write.fasta(sequences = write.align, names = names(write.align),file.out=locus.name.fa, nbchar = 1000000, as.string = T)
	## name of write.align (other than format, this is the same alignment as input.align)
	input.file      <- locus.name.fa                                                 ### Name of the temporary input fasta file
	output.file     <- paste0(input.file,"-tm")                                      ### Name of the temporary output trimal file
	trimAL.exe.path <- paste0(trimal.exe.dir,"/",trimal.exe.name)                    ### Full path to the trimal executable
	system(paste(trimAL.exe.path,"-in",input.file,"-out",output.file,"-automated1")) ### Executes trimAL
	system(paste("rm",input.file))                                                   ### Removes the fasta version of input.align
	
	## if trimAL alignment does not exist, delete the fasta version of input.align and print a warning
	if (file.exists(output.file) == F) { 
		system(paste0("rm ", input.file))
		print(paste0("deleted. Not enough overlapping data in alignment.") )
		return(Biostrings::DNAStringSet())
	} else {
		system(paste("mv",output.file,input.file))  ### renames the trimAL alignment filename
	}
	
	out.align <- Rsamtools::scanFa(Rsamtools::FaFile(input.file))       ### reads the trimAL alignment
	
	#Fixes any terrible NA names introduced by trimal
	new.names<-c()
	for (j in 1:length(names(out.align))){ 
		new.names[j] <- save.rownames[grep(pattern = names(out.align)[j], x = save.rownames)]
	}
	
	temp <- names(out.align)[is.na(names(out.align)) == T]
	if (length(temp) > 0){
		stop("there are NAs in the names")
	}
	names(out.align) <- new.names
	unlink(input.file) ### deletes the trimAL alignment (because its not in its final form yet)
	return(out.align)  ### the value of the function is out.align
}#end function




