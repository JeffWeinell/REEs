#' @title Convert Multiple Alignment Format to Fasta
#'
#' For each MAF alignment block, this function writes a fasta alignment file
#' 
#' @param maf Path to the input MAF file
#' @param outdir Path to directory where output fasta files should be written
#' @param samples Names of samples to keep. Default NULL means all samples are kept.
#' @param minwidth Number indicating the minimum alignment width to retain an alignment. Default 100.
#' @param minsamples Number indicating the minimum number of individuals in an alignment to keep the alignment. Default 4.
#' @param silent Logical indicating whether or not to omit printing to screen. Default TRUE.
#' @return NULL; alignments are written to outdir
#' @export maf2fasta
maf2fasta <- function(maf,outdir,samples=NULL,minwidth=100,minsamples=4,silent=F){
	data.lines  <- readLines(maf)
	lines.a     <- grep("^a$",data.lines)
	lines.s1    <- lines.a+1
	lines.slast <- c((lines.a-2)[-1],(length(data.lines)-1))
	# for(i in 1:length(lines.a)){
	for(i in 10001:length(lines.a)){
		lines.temp  <- data.lines[lines.s1[i]:lines.slast[i]]
		mat.temp    <- do.call(rbind,strsplit(lines.temp,split="\t"))[,-1,drop=F]
		idmat.temp  <- do.call(rbind,strsplit(mat.temp[,1,drop=F],split="\\."))
		mat.temp2   <- cbind(idmat.temp,mat.temp[,-1,drop=F])
		colnames(mat.temp2) <- c("SampleName","SampleContig","SourceStart","SourceEnd","SourceStrand","SourceSize","Sequence")
		if(!is.null(samples)){
			#match(samples,mat.temp2[,"SampleName"])
			if(any(mat.temp2[,"SampleName"] %in% samples)){
				mat.temp3 <- mat.temp2[which(mat.temp2[,"SampleName"] %in% samples),,drop=F]
			} else {
				next
			}
		} else {
			mat.temp3 <- mat.temp2
		}
		dna.temp <- Biostrings::DNAStringSet(mat.temp3[,"Sequence"])
		names(dna.temp) <- sprintf("Alignment%s_%s_%s",i,mat.temp3[,"SampleName"],mat.temp3[,"SampleContig"])
		if(any(width(dna.temp) >= minwidth)){
			dna.temp <- dna.temp[which(width(dna.temp) >= minwidth)]
		} else {
			next
		}
		if(length(dna.temp) < minsamples){
			next
		}
		Biostrings::writeXStringSet(dna.temp,file.path(outdir,paste0("alignment",i,".fa")))
		print(i)
	}
}

