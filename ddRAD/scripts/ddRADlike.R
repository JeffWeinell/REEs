#' In Silico ddRAD
#' 
#' Perform in silico ddRAD on a reference genome.
#' Note: not quite like doing ddRAD because instead of fragmenting and then size selecting, the method implemented here finds ranges on a contig between two different RE cut sites – but these ranges (potential targets) may have cut sites WITHIN them, without consequence.
#' 
#' @param input.seqs Either NULL (default but not preferred usage), or a character string indicating either a local filepath or URL path to the input sequences in fasta format. The URL character string is often the preferred method.
#' @param output.dir Character string indicating the directory where output files should be saved.
#' @param recognitionSeqs A character vector of length 2 containing the restriction enzyme recognition sequences to use.
#' @param lim.lengths Numerical vector of length 2, containing the minimum and maximum fragment lengths to keep.
#' @param save.tables Should the tables holding the locations of each recognition sequence on each contig be saved? Default is TRUE but these tables can be big.
#' @return A table containing the proposed loci and their contig coordinates. Output files (optionally) include two tables, one for each RE recognition sequence, holding the positions of recognition sequences in each genome contig.
#' @export proposeLoci.ddRADlike
proposeLoci.ddRADlike <- function(input.seqs,output.dir,recognitionSeqs,lim.lengths=c(500,1000),save.tables=T){
	if("DNAStringSet" %in% class(input.seqs)){
			fa0                 <- input.seqs
			delete.subject      <- F
	} else {
		if(file.exists(input.seqs)){
			fa0            <- Biostrings::readDNAStringSet(input.seqs)
			delete.subject <- F
		} else {
			### Define path that will be used to download the genome as a temporary file.
			subject.path <- tempfile()
			### Set time limit for downloading files to 1000 seconds (or longer if needed)
			options(timeout=1000)
			### Download genome to the temporary file
			utils::download.file(url=input.seqs, destfile=subject.path)
			fa0                  <- Biostrings::readDNAStringSet(subject.path)
			delete.subject       <- T
		}
	}
	headers              <- gsub(" .+","",names(fa0))
	names(fa0)           <- headers

	### Define restriction enzyme (RE) recognition sequences
	RecognitionSeqA.pattern <- recognitionSeqs[1]
	RecognitionSeqB.pattern <- recognitionSeqs[2]

	### Calculate lengths of RE recognition sequences
	RecognitionSeqA.length <- nchar(RecognitionSeqA.pattern)
	RecognitionSeqB.length <- nchar(RecognitionSeqB.pattern)

	### Find start positions of recognition sequence A on each contig
	RecognitionSeqA.start.pos        <- gregexpr(RecognitionSeqA.pattern, fa0)
	names(RecognitionSeqA.start.pos) <- names(fa0)
	RecognitionSeqA.start.pos        <- unlist(RecognitionSeqA.start.pos)
	
	### Drop entries with value "-1", which means the recognition site was not found on the contig.
	if(any(RecognitionSeqA.start.pos==-1)){
		RecognitionSeqA.start.pos <- RecognitionSeqA.start.pos[-which(RecognitionSeqA.start.pos== -1)]
	}
	
	### Calculate end positions of recognition sequence A on each contig
	RecognitionSeqA.end.pos   <- RecognitionSeqA.start.pos+RecognitionSeqA.length-1
	
	### Find start positions of recognition sequence B on each contig
	RecognitionSeqB.start.pos        <- gregexpr(RecognitionSeqB.pattern, fa0)
	names(RecognitionSeqB.start.pos) <- names(fa0)
	RecognitionSeqB.start.pos        <- unlist(RecognitionSeqB.start.pos)
	
	### Drop entries with value "-1", which means the recognition site was not found on the contig.
	if(any(RecognitionSeqB.start.pos==-1)){
		RecognitionSeqB.start.pos <- RecognitionSeqB.start.pos[-which(RecognitionSeqB.start.pos== -1)]
	}
	### End positions of Recognition Sequence B relative to query sequence (ith contig)
	RecognitionSeqB.end.pos   <- RecognitionSeqB.start.pos+RecognitionSeqB.length-1
	
	### Create a matrix of contigs, start positions, and end positions for RecognitionSeqA and RecognitionSeqB
	RecognitionSeqA.data <- cbind(names(RecognitionSeqA.start.pos),RecognitionSeqA.start.pos,RecognitionSeqA.end.pos)
	RecognitionSeqB.data <- cbind(names(RecognitionSeqB.start.pos),RecognitionSeqB.start.pos,RecognitionSeqB.end.pos)
	
	### Define column names for RecognitionSeqA.data and RecognitionSeqB.data
	colnames(RecognitionSeqA.data) <- c("contig.accession","RecognitionSeqA.start.pos","RecognitionSeqA.end.pos")
	colnames(RecognitionSeqB.data) <- c("contig.accession","RecognitionSeqB.start.pos","RecognitionSeqB.end.pos")
	
	### Update the Genbank accession numbers in column 1 of RecognitionSeqA.data and RecognitionSeqB.data.
	RecognitionSeqA.data[,"contig.accession"] <- gsub("\\.1.+$","\\.1",RecognitionSeqA.data[,"contig.accession"])
	RecognitionSeqB.data[,"contig.accession"] <- gsub("\\.1.+$","\\.1",RecognitionSeqB.data[,"contig.accession"])
	
	### Where to save RecognitionSeqA.data and RecognitionSeqB.data
	if(save.tables==T){
		outfile1 <- paste0(output.dir,RecognitionSeqA.pattern,"_RecognitionSeqATable.txt")
		outfile2 <- paste0(output.dir,RecognitionSeqB.pattern,"_RecognitionSeqBTable.txt")
	} else {
		outfile1 <- tempfile()
		outfile2 <- tempfile()
	}

	### Write RecognitionSeqA.data and RecognitionSeqB.data tables to files
	write.table(RecognitionSeqA.data,file=outfile1 ,sep="\t",append=F,col.names=T,row.names=F,quote=F)
	write.table(RecognitionSeqB.data,file=outfile2 ,sep="\t",append=F,col.names=T,row.names=F,quote=F)


	### Load in grep Hit Tables
#	RecognitionSeqATable <- read.table(file=paste0(directory,RecognitionSeqA.pattern,"_RecognitionSeqATable.txt"),header=T,sep="\t")
#	RecognitionSeqBTable <- read.table(file=paste0(directory,RecognitionSeqB.pattern,"_RecognitionSeqBTable.txt"),header=T,sep="\t")
	#rm(RecognitionSeqA.data)
	#rm(RecognitionSeqB.data)
	RecognitionSeqATable <- RecognitionSeqA.data
	RecognitionSeqBTable <- RecognitionSeqB.data
	# RecognitionSeqATable <- as.data.frame(data.table::fread(outfile1))
	# RecognitionSeqBTable <- as.data.frame(data.table::fread(outfile2))

#	output.name      <- paste0(directory,"proposedLoci_",RecognitionSeqA.pattern,"-",RecognitionSeqB.pattern,"_output.txt")
#	comp.output.name <- paste0(directory,"compProposedLoci_",RecognitionSeqA.pattern,"-",RecognitionSeqB.pattern,"_output.txt")
	if(save.tables==T){
		coordinates.output.name <- paste0(output.dir,"ProposedLoci_",RecognitionSeqA.pattern,"-",RecognitionSeqB.pattern,"_output.txt")
	} else {
		coordinates.output.name <- tempfile()
	}
	#accessionColumn <- 1
	#startPosColumn  <- 2
	#endPosColumn    <- 3

	min.length              <- lim.lengths[1]
	max.length              <- lim.lengths[2]

	# contigs <- unique(c(RecognitionSeqATable[,1], RecognitionSeqBTable[,1]))
	### Get the list of contigs that contain at least one of each of the recognition sites
	contigs <- intersect(RecognitionSeqATable[,1], RecognitionSeqBTable[,1])

	result <- list(); length(result) <- length(contigs)

	for(i in 1:length(contigs)){
		#if(i==1){
		#	headers=matrix(data=c("coordinates","ContigAccession","StartPosition","EndPosition","TargetLength"),nrow=1)
		#	write.table(headers,file=output.name,sep="\t",append=F,col.names=F,row.names=F,quote=F)
		#	write.table(headers,file=comp.output.name,sep="\t",append=F,col.names=F,row.names=F,quote=F)
		#}
		RecognitionSeqASet.temp <- which(RecognitionSeqATable[,1]==contigs[i])
		RecognitionSeqBSet.temp <- which(RecognitionSeqBTable[,1]==contigs[i])
		
		#RecognitionSeqA.startPositions <- as.numeric(RecognitionSeqATable[RecognitionSeqASet.temp,2])
		#RecognitionSeqA.endPositions   <- as.numeric(RecognitionSeqATable[RecognitionSeqASet.temp,3])
		#RecognitionSeqB.startPositions <- as.numeric(RecognitionSeqBTable[RecognitionSeqBSet.temp,2])
		#RecognitionSeqB.endPositions   <- as.numeric(RecognitionSeqBTable[RecognitionSeqBSet.temp,3])
	
	##### Added this code <------
		start.end.A             <- paste(RecognitionSeqATable[RecognitionSeqASet.temp,2],RecognitionSeqATable[RecognitionSeqASet.temp,3],sep="_")
		start.end.B             <- paste(RecognitionSeqBTable[RecognitionSeqBSet.temp,2],RecognitionSeqBTable[RecognitionSeqBSet.temp,3],sep="_")
		startA.endA.startB.endB <- mat.strsplit(gsub("_"," ",xprod.combn(start.end.A,start.end.B)))
		mode(startA.endA.startB.endB) <- "numeric"
		loci.start      <- apply(X=startA.endA.startB.endB,MARGIN=1,FUN=function(row){min(row[c(1,3)])})
		loci.end        <- apply(X=startA.endA.startB.endB,MARGIN=1,FUN=function(row){max(row[c(2,4)])})
		lengths         <- loci.end-loci.start+1
		#loci.mat        <- cbind(rep(contigs[1],length(keep1)),loci.start,loci.end,lengths)
		#loci.ranges     <- cbind(loci.start,loci.end)
		keep1           <- which(lengths > min.length & lengths < max.length)
		if(length(keep1)==0){
			next
		}
		ContigAccession <- rep(contigs[i],length(keep1))
		StartPosition   <- loci.start[keep1]
		EndPosition     <- loci.end[keep1]
		TargetLength    <- lengths[keep1]
		take.rv.Comp    <- (startA.endA.startB.endB[,1] > startA.endA.startB.endB[,3])[keep1]
		coordinates     <- paste0(ContigAccession,":",StartPosition,"-",EndPosition)
		coordinates[take.rv.Comp] <- paste0(ContigAccession,":c",EndPosition,"-",StartPosition)[take.rv.Comp]
		targets.mat.temp <- cbind(coordinates,ContigAccession,StartPosition,EndPosition,TargetLength)
		result[[i]] <- targets.mat.temp
	
	##### ----> End of new block of code.
	
	#	distances      <- matrix(0, nrow=length(RecognitionSeqASet.temp), ncol=length(RecognitionSeqBSet.temp))
	#	for(j in 1:length(RecognitionSeqASet.temp)){
	#		starts.temp    <- rep(RecognitionSeqA.startPositions[j],length(RecognitionSeqBSet.temp))
	#		#lengths.temp  <- RecognitionSeqB.endPositions-starts.temp
	#		lengths.temp   <- (RecognitionSeqB.endPositions-starts.temp)+1
	#		distances[j,]  <- lengths.temp
	#	}
	#	
	#	comp.distances <- matrix(0, nrow=length(RecognitionSeqBSet.temp), ncol=length(RecognitionSeqASet.temp))
	#	for(z in 1:length(RecognitionSeqBSet.temp)){
	#		comp.starts.temp      <- rep(RecognitionSeqB.startPositions[z], length(RecognitionSeqASet.temp))
	#		#comp.lengths.temp     <- RecognitionSeqA.endPositions-comp.starts.temp
	#		comp.lengths.temp     <- (RecognitionSeqA.endPositions-comp.starts.temp)+1
	#		comp.distances[z,]    <- comp.lengths.temp
	#	}
	#	
	#	matches.proposed        <- which(distances > min.length & distances < max.length,arr.ind=TRUE)
	#	comp.matches.proposed   <- which(comp.distances > min.length & comp.distances < max.length,arr.ind=TRUE)
	#	if(nrow(matches.proposed)!=0){
	#		matches.mat           <- matrix(data=0,nrow=nrow(matches.proposed),ncol=5)
	#		colnames(matches.mat) <- c("coordinates","ContigAccession","StartPosition","EndPosition","TargetLength")
	#		for(k in 1:nrow(matches.proposed)){
	#			start.position.temp <- RecognitionSeqA.startPositions[matches.proposed[k,1]]
	#			end.position.temp   <- RecognitionSeqB.endPositions[matches.proposed[k,2]]
	##			accession.temp      <- contigs[i]
	#			coordinates.temp    <- paste0(contigs[i],":",start.position.temp,"-",end.position.temp)
	#			distance.temp       <- distances[matches.proposed[k,1],matches.proposed[k,2]]
	#			info.temp           <- matrix(data=c(coordinates.temp,contigs[i],start.position.temp,end.position.temp,distance.temp),nrow=1)
	#			matches.mat[k,]     <- info.temp
	##			write.table(info.temp,file=output.name,sep="\t",append=T,col.names=F,row.names=F,quote=F)
	#		}
	#	}
	#	
	#	if(nrow(comp.matches.proposed)!=0){
	#		comp.matches.mat           <- matrix(data=0,nrow=nrow(comp.matches.proposed),ncol=5)
	#		colnames(comp.matches.mat) <- c("coordinates","ContigAccession","StartPosition","EndPosition","TargetLength")
	#		for(x in 1:nrow(comp.matches.proposed)){
	#			comp.start.position.temp <- RecognitionSeqB.startPositions[comp.matches.proposed[x,1]]
	#			comp.end.position.temp   <- RecognitionSeqA.endPositions[comp.matches.proposed[x,2]]
	##			comp.accession.temp      <- contigs[i]
	#			comp.coordinates.temp    <- paste(contigs[i],":c",comp.end.position.temp,"-",comp.start.position.temp,sep="")
	#			comp.distance.temp       <- comp.distances[comp.matches.proposed[x,1],comp.matches.proposed[x,2]]
	#			comp.info.temp           <- matrix(data=c(comp.coordinates.temp,contigs[i],comp.start.position.temp,comp.end.position.temp,comp.distance.temp),nrow=1)
	#			comp.matches.mat[x,]     <- c(comp.info.temp)
	##			write.table(comp.info.temp,file=comp.output.name,sep="\t",append=T,col.names=F,row.names=F,quote=F)
	#		}
	#	}	
	}
	result <- do.call(rbind,result)
	write.table(result,file=coordinates.output.name,sep="\t",append=F,col.names=T,row.names=F,quote=F)
	result
}


#	PstI   <- "CTGCAG"
#	HpaII  <- "CCGG"
#	SbfI   <- "CCTGCAGG"
#	EcoRI  <- "GAATTC"
#	MluCI  <- "AATT"
#	NlaIII <- "CATG"
#
#	RecognitionSeqA.pattern <- SbfI
#	RecognitionSeqB.pattern <- EcoRI

