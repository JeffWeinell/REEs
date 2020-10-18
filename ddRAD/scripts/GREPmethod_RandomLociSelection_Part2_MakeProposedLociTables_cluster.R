directory <- "~/workingDirName"

PstI   <- "CTGCAG"
HpaII  <- "CCGG"
SbfI   <- "CCTGCAGG"
EcoRI  <- "GAATTC"
MluCI  <- "AATT"
NlaIII <- "CATG"

RecognitionSeqA.pattern <- SbfI
RecognitionSeqB.pattern <- EcoRI

### Load in grep Hit Tables
RecognitionSeqATable <- read.table(file=paste0(directory,RecognitionSeqA.pattern,"_RecognitionSeqATable.txt"),header=T,sep="\t")
RecognitionSeqBTable <- read.table(file=paste0(directory,RecognitionSeqB.pattern,"_RecognitionSeqBTable.txt"),header=T,sep="\t")

output.name      <- paste0(directory,"proposedLoci_",RecognitionSeqA.pattern,"-",RecognitionSeqB.pattern,"_output.txt")
comp.output.name <- paste0(directory,"compProposedLoci_",RecognitionSeqA.pattern,"-",RecognitionSeqB.pattern,"_output.txt")

accessionColumn <- 1
startPosColumn  <- 2
endPosColumn    <- 3

Contigs <- unique(c(as.character(RecognitionSeqBTable[,accessionColumn]), as.character(RecognitionSeqATable[,accessionColumn])))

for(i in 1:length(Contigs)){
	if(i==1){
		headers=matrix(data=c("coordinates","ContigAccession","StartPosition","EndPosition","TargetLength"),nrow=1)
		write.table(headers,file=output.name,sep="\t",append=F,col.names=F,row.names=F,quote=F)
		write.table(headers,file=comp.output.name,sep="\t",append=F,col.names=F,row.names=F,quote=F)
	}
	RecognitionSeqASet.temp <- which(as.character(RecognitionSeqATable[,accessionColumn])==Contigs[i])
	RecognitionSeqBSet.temp <- which(as.character(RecognitionSeqBTable[,accessionColumn])==Contigs[i])
	if(length(RecognitionSeqASet.temp)==0){
    next
   }
	if(length(RecognitionSeqBSet.temp)==0){
    next
   }
	
	RecognitionSeqA.startPositions <- as.numeric(RecognitionSeqATable[RecognitionSeqASet.temp,startPosColumn])
	RecognitionSeqA.endPositions   <- as.numeric(RecognitionSeqATable[RecognitionSeqASet.temp,endPosColumn])
	RecognitionSeqB.startPositions <- as.numeric(RecognitionSeqBTable[RecognitionSeqBSet.temp,startPosColumn])
	RecognitionSeqB.endPositions   <- as.numeric(RecognitionSeqBTable[RecognitionSeqBSet.temp,endPosColumn])
	
	distances      <- matrix(0, nrow=length(RecognitionSeqASet.temp), ncol=length(RecognitionSeqBSet.temp))
	comp.distances <- matrix(0, nrow=length(RecognitionSeqBSet.temp), ncol=length(RecognitionSeqASet.temp))
	
	for(j in 1:length(RecognitionSeqASet.temp)){
		starts.temp   <- rep(RecognitionSeqA.startPositions[j],length(RecognitionSeqBSet.temp))
		lengths.temp  <- RecognitionSeqB.endPositions-starts.temp
		distances[j,] <- lengths.temp
	}
	
	for(z in 1:length(RecognitionSeqBSet.temp)){
		comp.starts.temp      <- rep(RecognitionSeqB.startPositions[z], length(RecognitionSeqASet.temp))
		comp.lengths.temp     <- RecognitionSeqA.endPositions-comp.starts.temp
		comp.distances[z,]    <- comp.lengths.temp
	}
	
	min.length              <- 500
	max.length              <- 1000
	matches.proposed        <- which(distances > min.length & distances < max.length,arr.ind=TRUE)
	comp.matches.proposed   <- which(comp.distances > min.length & comp.distances < max.length,arr.ind=TRUE)
	
	if(nrow(matches.proposed)!=0){
		for(k in 1:nrow(matches.proposed)){
			start.position.temp <- RecognitionSeqA.startPositions[matches.proposed[k,1]]
			end.position.temp   <- RecognitionSeqB.endPositions[matches.proposed[k,2]]
			accession.temp      <- Contigs[i]
			coordinates.temp    <- paste(Contigs[i],":",start.position.temp,"-",end.position.temp,sep="")
			distance.temp       <- distances[matches.proposed[k,1],matches.proposed[k,2]]
			info.temp           <- matrix(data=c(coordinates.temp,accession.temp,start.position.temp,end.position.temp,distance.temp),nrow=1)
			write.table(info.temp,file=output.name,sep="\t",append=T,col.names=F,row.names=F,quote=F)
		}
	}
	if(nrow(comp.matches.proposed)!=0){
		for(x in 1:nrow(comp.matches.proposed)){
			comp.start.position.temp <- RecognitionSeqB.startPositions[comp.matches.proposed[x,1]]
			comp.end.position.temp   <- RecognitionSeqA.endPositions[comp.matches.proposed[x,2]]
			comp.accession.temp      <- Contigs[i]
			comp.coordinates.temp    <- paste(Contigs[i],":c",comp.end.position.temp,"-",comp.start.position.temp,sep="")
			comp.distance.temp       <- comp.distances[comp.matches.proposed[x,1],comp.matches.proposed[x,2]]
			comp.info.temp           <- matrix(data=c(comp.coordinates.temp,comp.accession.temp,comp.start.position.temp,comp.end.position.temp,comp.distance.temp),nrow=1)
			write.table(comp.info.temp,file=comp.output.name,sep="\t",append=T,col.names=F,row.names=F,quote=F)
		}
	}	
}
