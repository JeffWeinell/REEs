#' @title translate CDS regions
#' 
#' Used to translate the output of get.seqs.from.gff. Translates reading frame 1 of forward and reverse complement sequences of input DNAStringSet and returns whichever has fewer stop codons.
#' The output of this function can be the query input for blast function with method "tblastn"
#' 
#' @param input.seqs DNAStringSet produced by get.seqs.from.gff
#' @param drop.if.stop Should AA sequences with internal stop codons be dropped? Default is TRUE.
#' @return AAStringSet containing translated sequences of input.seqs
#' @export translate.exome
translate.exome <- function(input.seqs,drop.if.stop=T){
	### Obtain Forward and Reverse Complement reading frames (three each)
	F1          <- input.seqs
	F2          <- subseq(F1,start=2)
	F3          <- subseq(F1,start=3)
	RC1         <- Biostrings::reverseComplement(input.seqs)
	RC2         <- subseq(RC1,start=2)
	RC3         <- subseq(RC1,start=3)
	### Append frame info to sequence names
	names(F1)   <- paste0(names(F1),"_F1")
	names(F2)   <- paste0(names(F2),"_F2")
	names(F3)   <- paste0(names(F3),"_F3")
	names(RC1)  <- paste0(names(RC1),"_RC1")
	names(RC2)  <- paste0(names(RC2),"_RC2")
	names(RC3)  <- paste0(names(RC3),"_RC3")
	### Calculate subsequence widths that would be divisible by three
	widths.F1   <- floor(width(F1)/3)*3
	widths.F2   <- floor(width(F2)/3)*3
	widths.F3   <- floor(width(F3)/3)*3
	widths.RC1  <- floor(width(RC1)/3)*3
	widths.RC2  <- floor(width(RC2)/3)*3
	widths.RC3  <- floor(width(RC3)/3)*3
	### Update F1, F2, F3, RC1, RC2, RC3 such that each is divisible by three
	F1          <- subseq(F1,start=1,end=widths.F1)
	F2          <- subseq(F2,start=1,end=widths.F2)
	F3          <- subseq(F3,start=1,end=widths.F3)
	RC1         <- subseq(RC1,start=1,end=widths.RC1)
	RC2         <- subseq(RC2,start=1,end=widths.RC2)
	RC3         <- subseq(RC3,start=1,end=widths.RC3)
	### Create a list containing the translated reading frame of each DNA sequence in input.seqs with the fewest number of stop codons. Take the first among ties.
#	result.list <- lapply(X=c(1:length(input.seqs)),FUN=function(x){temp.seqs = suppressWarnings(Biostrings::translate(c(F1[x],F2[x],F3[x],RC1[x],RC2[x],RC3[x]))); counts.temp = REEs::countChars(as.character(temp.seqs),pattern="\\*"); temp.seqs[which(counts.temp==min(counts.temp))]})
	result.list   <- lapply(X=c(1:length(input.seqs)),FUN=function(x){temp.seqs = suppressWarnings(Biostrings::translate(c(F1[x],F2[x],F3[x],RC1[x],RC2[x],RC3[x]))); temp.seqs})
	### Collapse the list of AAStringSets into a single AAStringSet
	result        <- REEs::collapse.DNAStringSet(result.list)
	names(result) <- gsub("_$","",names(result))
	if(drop.if.stop){
		result    <- result[-grep("\\*.",test)]
	}
	### Return result
	result
}





