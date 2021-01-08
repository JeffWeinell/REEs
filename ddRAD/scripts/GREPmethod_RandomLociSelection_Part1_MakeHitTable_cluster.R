directory="~/WorkingDirName"

files <- list.files(path="~/Thermophis_Genome_Contigs/")

RecognitionSeqA.pattern <- "CTGCAG"    ## PstI recognition site
RecognitionSeqB.pattern <- "CCGG"      ## HpaII recognition site

RecognitionSeqA.length <- nchar(RecognitionSeqA.pattern) ### nucleotide length of Recognition Sequence A
RecognitionSeqB.length <- nchar(RecognitionSeqB.pattern) ### nucleotide length of Recognition Sequence B

### Name of output files
outfile1 <- paste0(directory,RecognitionSeqA.pattern,"_RecognitionSeqATable.txt")
outfile2 <- paste0(directory,RecognitionSeqB.pattern,"_RecognitionSeqBTable.txt")


for(i in 1:length(files)){
	if(i==1){
		headers=matrix(data=c("accession","start","end"),nrow=1)
		write.table(headers,file=outfile1 ,sep="\t",append=F,col.names=F,row.names=F,quote=F)
		write.table(headers,file=outfile2 ,sep="\t",append=F,col.names=F,row.names=F,quote=F)
	}
	accession     <- gsub(".fasta","",files[i]) ### accession ID of ith contig 
	filename.temp <- paste0(directory,files[i]) ### full filepath to ith contig
	### Read in ith contig as a table
	sequences     <- read.table(file=filename.temp,header=T,sep="\t")
	### Collapse bases of ith contig into a character string.
	sequence.temp <- paste(as.character(unlist(sequences)),collapse="")

	RecognitionSeqA.start.pos <- unlist(gregexpr(RecognitionSeqA.pattern, sequence.temp))		## start positions of Recognition Sequence A relative to query sequence (ith contig)
	RecognitionSeqA.end.pos   <- RecognitionSeqA.start.pos+RecognitionSeqA.length-1			## end positions of Recognition Sequence A relative to query sequence (ith contig)
	RecognitionSeqB.start.pos <- unlist(gregexpr(RecognitionSeqB.pattern, sequence.temp))           ## start positions of Recognition Sequence B relative to query sequence (ith contig)
	RecognitionSeqB.end.pos   <- RecognitionSeqB.start.pos+RecognitionSeqB.length-1                 ## end positions of Recognition Sequence B relative to query sequence (ith contig)
	
	RecognitionSeqA.data <- cbind(accession,RecognitionSeqA.start.pos,RecognitionSeqA.end.pos)      ## a vector containing the info for Recognition Sequence A hits relative to the query sequence
	RecognitionSeqB.data <- cbind(accession,RecognitionSeqB.start.pos,RecognitionSeqB.end.pos)      

	write.table(RecognitionSeqA.data,file=outfile1 ,sep="\t",append=T,col.names=F,row.names=F,quote=F)
	write.table(RecognitionSeqB.data,file=outfile2,sep="\t",append=T,col.names=F,row.names=F,quote=F)
}
