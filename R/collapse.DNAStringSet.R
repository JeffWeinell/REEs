#' @title Collapse DNAStringSets
#' 
#' Collapse a list of DNAStringSet objects into a single DNAStringSet object.
#' 
#' @param input.seqs A list of DNAStringSet objects.
#' @param use.seqnames Should the names of the input DNAStringSets be incorporated into the seqnames of the collapsed DNAStringSet. Default T.
#' @param use.set.names Should the names of the input DNAStringSets be incorporated into the seqnames of the collapsed DNAStringSet. Default T.
#' @return A single DNAStringSet object containing the sequences from all of the DNAStringSet objects in input.seqs. Sequence names are determined by use.seqnames and use.set.names parameters.
#' @export collapse.DNAStringSet
collapse.DNAStringSet <- function(input.seqs,use.seqnames=T,use.set.names=T){
	names.original    <- lapply(input.seqs,FUN=names)
	seqnames.original <- unlist(names.original)
	seqs.temp         <- lapply(input.seqs,function(x){names(x)=NULL;x})
	names(seqs.temp)  <- NULL
	seqs.temp2        <- do.call(c,seqs.temp)
	if(use.seqnames & use.set.names){
		seqnames.out  <- paste0(seqnames.original,"_",names(seqnames.original))
	}
	if(use.seqnames & !use.set.names){
		seqnames.out  <- seqnames.original
	}
	if(!use.seqnames & use.set.names){
		seqnames.out  <- names(seqnames.original)
	}
	if(!use.seqnames & !use.set.names){
		seqnames.out  <- NULL
	}
	names(seqs.temp2) <- seqnames.out
	seqs.temp2
}