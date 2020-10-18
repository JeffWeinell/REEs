directory="/panfs/pfs.local/home/j926w878/scratch/"
files <- list.files(path="/panfs/pfs.local/home/j926w878/scratch/Thermophis_Genome_Contigs/")

RecognitionSeqA.pattern <- "CTGCAG"    ## PstI recognition site
RecognitionSeqB.pattern <- "CCGG"      ## HpaII recognition site

RecognitionSeqA.length <- nchar(RecognitionSeqA.pattern) ### nucleotide length of Recognition Sequence A
RecognitionSeqB.length <- nchar(RecognitionSeqB.pattern) ### nucleotide length of Recognition Sequence B

outfile1 <- paste(directory,RecognitionSeqA.pattern,"_RecognitionSeqATable.txt",sep="")
outfile2 <- paste(directory,RecognitionSeqB.pattern,"_RecognitionSeqBTable.txt",sep="")


for(i in 1:length(files)){
	if(i==1){
		headers=matrix(data=c("accession","start","end"),nrow=1)
		write.table(headers,file=outfile1 ,sep="\t",append=F,col.names=F,row.names=F,quote=F)
		write.table(headers,file=outfile2 ,sep="\t",append=F,col.names=F,row.names=F,quote=F)
	}
	accession <- gsub(".fasta","",files[i])
	filename.temp <- paste(directory,files[i],sep="")
	sequences <- read.table(file=filename.temp,header=T,sep="\t")
	sequence.temp <- paste(as.character(unlist(sequences)),collapse="")

	RecognitionSeqA.start.pos <- unlist(gregexpr(RecognitionSeqA.pattern, sequence.temp))			## start positions of Recognition Sequence A relative to query sequence
	RecognitionSeqA.end.pos <- RecognitionSeqA.start.pos+RecognitionSeqA.length-1					## end positions of Recognition Sequence A relative to query sequence
	RecognitionSeqB.start.pos <- unlist(gregexpr(RecognitionSeqB.pattern, sequence.temp))           ## start positions of Recognition Sequence B relative to query sequence
	RecognitionSeqB.end.pos <- RecognitionSeqB.start.pos+RecognitionSeqB.length-1                   ## end positions of Recognition Sequence B relative to query sequence
	
	RecognitionSeqA.data <- cbind(accession,RecognitionSeqA.start.pos,RecognitionSeqA.end.pos)      ## a vector containing the info for Recognition Sequence A hits relative to the query sequence
	RecognitionSeqB.data <- cbind(accession,RecognitionSeqB.start.pos,RecognitionSeqB.end.pos)      

	write.table(RecognitionSeqA.data,file=outfile1 ,sep="\t",append=T,col.names=F,row.names=F,quote=F)
	write.table(RecognitionSeqB.data,file=outfile2,sep="\t",append=T,col.names=F,row.names=F,quote=F)
}
