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
proposeLoci.ddRADlike <- function(input.seqs,output.dir,recognitionSeqs,lim.lengths=c(900,1000),save.tables=T){
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

	if(save.tables==T){
		coordinates.output.name <- paste0(output.dir,"ProposedLoci_",RecognitionSeqA.pattern,"-",RecognitionSeqB.pattern,"_output.txt")
	} else {
		coordinates.output.name <- tempfile()
	}
	### Get the list of contigs that contain at least one of each of the recognition sites
	contigs <- intersect(RecognitionSeqATable[,1], RecognitionSeqBTable[,1])
	startA  <- apply(RecognitionSeqATable[,2:3],1,min)
	endA    <- apply(RecognitionSeqATable[,2:3],1,max)
	startB  <- apply(RecognitionSeqBTable[,2:3],1,min)
	endB    <- apply(RecognitionSeqBTable[,2:3],1,max)

	startA.list <- lapply(contigs,FUN=function(x){startA[which(RecognitionSeqATable[,"contig.accession"]==x)]})
	endA.list   <- lapply(contigs,FUN=function(x){endA[which(RecognitionSeqATable[,"contig.accession"]==x)]})
	startB.list <- lapply(contigs,FUN=function(x){startB[which(RecognitionSeqBTable[,"contig.accession"]==x)]})
	endB.list   <- lapply(contigs,FUN=function(x){endB[which(RecognitionSeqBTable[,"contig.accession"]==x)]})
	names(startA.list) <- contigs
	names(endA.list)   <- contigs
	names(startB.list) <- contigs
	names(endB.list)   <- contigs

	#AB.keep <- mapply(FUN=function(A,B){},A=startA.list,B=endB.list)
	#BA.keep <- mapply(FUN=function(A,B){},A=startB.list,B=endA.list)

	function1 <- function(A,min.len=900,max.len=1000){c(sapply(A,FUN=function(x){(x+min.len-1):(x+max.len-1)}))}
	function2 <- function(A,min.len=900,max.len=1000){c(lapply(A,FUN=function(x){(x+min.len-1):(x+max.len-1)}))}

	#any(B %in% function1(A))
	keep.contigs.AB <- mapply(function(A,B){any(B %in% function1(A))},A=startA.list,B=endB.list)
	keep.contigs.BA <- mapply(function(A,B){any(B %in% function1(A))},A=startB.list,B=endA.list)

	### For contigs with recognition site B that is 900-1000 downstream of a recognition site A, return which recognition site B is 900-1000 downstream of a recognition site A.
	endB.list.keep     <- mapply(function(A,B){which(B %in% function1(A))},A=startA.list[keep.contigs.AB],B=endB.list[keep.contigs.AB])
	### Keep recognition site B positions that are 900-1000 downstream of a recognition site A.
	endB.list.filtered <- mapply(function(A,B){A[B]},A=endB.list[keep.contigs.AB],B=endB.list.keep)

	### For contigs with recognition site A that is 900-1000 downstream of a recognition site B, return which recognition site A is 900-1000 downstream of a recognition site B.
	endA.list.keep     <- mapply(function(A,B){which(B %in% function1(A))},A=startB.list[keep.contigs.BA],B=endA.list[keep.contigs.BA])
	### Keep recognition site B positions that are 900-1000 downstream of a recognition site A.
	endA.list.filtered <- mapply(function(A,B){A[B]},A=endA.list[keep.contigs.BA],B=endA.list.keep)

	### For contigs with recognition site B that is 900-1000 downstream of a recognition site A, return which recognition site A is 900-1000 upstream of a recognition site B.
	startA.list.keep <- mapply(function(A,B){which(lapply(function2(A),FUN=function(x){any(x %in% B)})==T)},A=startA.list[keep.contigs.AB],B=endB.list[keep.contigs.AB])
	### Keep recognition site A positions that are 900-1000 upstream of a recognition site B.
	startA.list.filtered <- mapply(function(A,B){A[B]},A=startA.list[keep.contigs.AB],B=startA.list.keep)

	### For contigs with recognition site A that is 900-1000 downstream of a recognition site B, return which recognition site B is 900-1000 upstream of a recognition site A.
	startB.list.keep <- mapply(function(A,B){which(lapply(function2(A),FUN=function(x){any(x %in% B)})==T)},A=startB.list[keep.contigs.BA],B=endA.list[keep.contigs.BA])
	### Keep recognition site B positions that are 900-1000 upstream of a recognition site A.
	startB.list.filtered <- mapply(function(A,B){A[B]},A=startB.list[keep.contigs.BA],B=startB.list.keep)

#	keep.contigs    <- mapply(FUN=any,keep.contigs.AB,keep.contigs.BA)
	contigs.AB <- contigs[keep.contigs.AB]
	contigs.BA <- contigs[keep.contigs.BA]
	
	


#	RecognitionSeqASet <- lapply(X=contigs,FUN=function(x){which(RecognitionSeqATable[,1]==x)})
#	RecognitionSeqBSet <- lapply(X=contigs,FUN=function(x){which(RecognitionSeqBTable[,1]==x)})
#	start.end.A        <- lapply(X=RecognitionSeqASet,FUN=function(x){paste(RecognitionSeqATable[x,2],RecognitionSeqATable[x,3],sep="_")})
#	start.end.B        <- lapply(X=RecognitionSeqBSet,FUN=function(x){paste(RecognitionSeqBTable[x,2],RecognitionSeqBTable[x,3],sep="_")})
#	which()
#	any(start.end.B[[1]] %in% start.end.A[[1]]-1000)

	#startA.endA.startB.endB.temp  <- mapply(function(x,y){xprod.combn(x,y)},start.end.A,start.end.B)
	#startA.endA.startB.endB.temp  <- unlist(startA.endA.startB.endB.temp)
	#startA.endA.startB.endB.temp2 <- gsub("_"," ",startA.endA.startB.endB.temp)
	#startA.endA.startB.endB       <- mat.strsplit(startA.endA.startB.endB.temp2)
	#mode(startA.endA.startB.endB) <- "numeric"
	#hits.per.accessionA  <- sapply(start.end.A,FUN=length)
	#hits.per.accessionB  <- sapply(start.end.B,FUN=length)
	#combn.per.accession  <- hits.per.accessionA*hits.per.accessionB
	#combn.per.accession <- mapply(FUN=function(),start.end.A,start.end.B)
	temp.function <- function(A,B,min.len,max.len){
		#res.mat       <- mat.strsplit(gsub("_"," ",xprod.combn(A,B)))
		res.mat       <- mat.strsplit(xprod.combn(A,B))
		mode(res.mat) <- "numeric"
		#res.start     <- apply(X=res.mat,MARGIN=1,FUN=min)
		#res.end       <- apply(X=res.mat,MARGIN=1,FUN=max)
		#res.lengths   <- (res.end-res.start)+1
		res.lengths   <- res.mat[,2]-res.mat[,1]+1
		keep1         <- intersect(which(res.lengths >= min.len), which(res.lengths <= max.len))
		if(length(keep1)==0){
			res.mat <- NULL
		} else {
			res.mat <- res.mat[keep1,,drop=F]
		}
		res.mat
	}

	### Get a table of ranges for loci in which recognition site A is 900-1000nt upstream of recognition site B
	i.start.AB     <- 1
	i.end.AB       <- length(contigs.AB)
	startA.endB <- NULL
	pb = txtProgressBar(min = i.start.AB, max = i.end.AB, initial = 0)
	for(i in i.start.AB:i.end.AB){
		startA.endB.i <- temp.function(A=startA.list.filtered[[i]],B=endB.list.filtered[[i]],min.len=lim.lengths[1],max.len=lim.lengths[2])
		#startA.endA.startB.endB.i <- temp.function(x=start.end.A[[i]],y=start.end.B[[i]],min=1000,max=1000)
		if(is.null(startA.endB.i)){
			setTxtProgressBar(pb,i)
			next
		} else {
			rownames(startA.endB.i) <- paste0(rep(contigs.AB[i],nrow(startA.endB.i)),1:nrow(startA.endB.i))
			#if(i==i.start){
			#	startA.endA.startB.endB <- startA.endA.startB.endB.i
			#} else {
				#startA.endA.startB.endB <- do.call(rbind,startA.endA.startB.endB.i)
			startA.endB <- rbind(startA.endB,startA.endB.i)
			#}
		}
		setTxtProgressBar(pb,i)
	}
	ContigAccession.AB <- gsub("\\.1.+$","\\.1",rownames(startA.endB))
	coordinates.AB     <- paste0(ContigAccession.AB,":",startA.endB[,1],"-",startA.endB[,2])
	lengths.AB         <- as.numeric(startA.endB[,2])-as.numeric(startA.endB[,1])
	result.AB          <- cbind(coordinates.AB,ContigAccession.AB,startA.endB,lengths.AB)

	### Get a table of ranges for loci in which recognition site B is 900-1000nt upstream of recognition site A
	i.start.BA  <- 1
	i.end.BA    <- length(contigs.BA)
	startB.endA <- NULL
	pb = txtProgressBar(min = i.start.BA, max = i.end.BA, initial = 0)
	for(i in i.start.BA:i.end.BA){
		startB.endA.i <- temp.function(A=startB.list.filtered[[i]],B=endA.list.filtered[[i]],min.len=lim.lengths[1],max.len=lim.lengths[2])
		#startA.endA.startB.endB.i <- temp.function(x=start.end.A[[i]],y=start.end.B[[i]],min=1000,max=1000)
		if(is.null(startB.endA.i)){
			setTxtProgressBar(pb,i)
			next
		} else {
			rownames(startB.endA.i) <- paste0(rep(contigs.BA[i],nrow(startB.endA.i)),1:nrow(startB.endA.i))
			#if(i==i.start){
			#	startA.endA.startB.endB <- startA.endA.startB.endB.i
			#} else {
				#startA.endA.startB.endB <- do.call(rbind,startA.endA.startB.endB.i)
			startB.endA <- rbind(startB.endA,startB.endA.i)
			#}
		}
		setTxtProgressBar(pb,i)
	}
	ContigAccession.BA <- gsub("\\.1.+$","\\.1",rownames(startB.endA))
	coordinates.BA     <- paste0(ContigAccession.BA,":c",startB.endA[,2],"-",startB.endA[,1])
	lengths.BA         <- as.numeric(startB.endA[,2])-as.numeric(startB.endA[,1])
	result.BA          <- cbind(coordinates.BA,ContigAccession.BA,startB.endA,lengths.BA)

#	#startA.endA.startB.endB <- mapply(FUN=function(x,y){res=mat.strsplit(gsub("_"," ",xprod.combn(x,y)));mode(res)="numeric";res},start.end.A,start.end.B)
#	#startA.endA.startB.endB <- mapply(FUN=temp.function,x=start.end.A[[1:5]],y=start.end.B[[1:5]],min=900,max=1000)
#	startA.endA.startB.endB <- list(); length(startA.endA.startB.endB) <- length(contigs)
#	i.start <- 1
#	i.end   <- length(contigs)
#	startA.endA.startB.endB <- NULL
#	pb = txtProgressBar(min = i.start, max = i.end, initial = 0)
#	for(i in i.start:i.end){
#		startA.endA.startB.endB.i <- temp.function(A=start.end.A[[i]],B=start.end.B[[i]],min.len=lim.lengths[1],max.len=lim.lengths[2])
#		#startA.endA.startB.endB.i <- temp.function(x=start.end.A[[i]],y=start.end.B[[i]],min=1000,max=1000)
#		if(is.null(startA.endA.startB.endB.i)){
#			setTxtProgressBar(pb,i)
#			next
#		} else {
#			rownames(startA.endA.startB.endB.i) <- paste0(rep(contigs[i],nrow(startA.endA.startB.endB.i)),1:nrow(startA.endA.startB.endB.i))
#			#if(i==i.start){
#			#	startA.endA.startB.endB <- startA.endA.startB.endB.i
#			#} else {
#				#startA.endA.startB.endB <- do.call(rbind,startA.endA.startB.endB.i)
#			startA.endA.startB.endB <- rbind(startA.endA.startB.endB,startA.endA.startB.endB.i)
#			#}
#		}
#		setTxtProgressBar(pb,i)
#	}
#	names(startA.endA.startB.endB) <- contigs
#	startA.endA.startB.endB        <- do.call(rbind,startA.endA.startB.endB)
#	ContigAccession                <- gsub("\\.1.+$","\\.1",rownames(startA.endA.startB.endB))

	###
	#result <- list(); length(result) <- length(contigs)
	#### Set up a progress bar for the for loop
	#pb = txtProgressBar(min = 0, max = length(contigs), initial = 0)
	#### Run the for loop
	#for(i in 1:length(contigs)){
	#	RecognitionSeqASet.temp <- which(RecognitionSeqATable[,1]==contigs[i])
	#	RecognitionSeqBSet.temp <- which(RecognitionSeqBTable[,1]==contigs[i])
	#	start.end.A             <- paste(RecognitionSeqATable[RecognitionSeqASet.temp,2],RecognitionSeqATable[RecognitionSeqASet.temp,3],sep="_")
	#	start.end.B             <- paste(RecognitionSeqBTable[RecognitionSeqBSet.temp,2],RecognitionSeqBTable[RecognitionSeqBSet.temp,3],sep="_")
	#	startA.endA.startB.endB <- mat.strsplit(gsub("_"," ",xprod.combn(start.end.A,start.end.B)))
	#	mode(startA.endA.startB.endB) <- "numeric"
	#	loci.start      <- apply(X=startA.endA.startB.endB,MARGIN=1,FUN=function(row){min(row[c(1,3)])})
	#	loci.end        <- apply(X=startA.endA.startB.endB,MARGIN=1,FUN=function(row){max(row[c(2,4)])})
	#	lengths         <- loci.end-loci.start+1
	#	keep1           <- which(lengths >= lim.lengths[1] & lengths <= lim.lengths[2])
	#	if(length(keep1)==0){
	#		next
	#	}
	#	ContigAccession  <- rep(contigs[i],length(keep1))
	#	StartPosition    <- loci.start[keep1]
	#	EndPosition      <- loci.end[keep1]
	#	TargetLength     <- lengths[keep1]
	#	take.rv.Comp     <- (startA.endA.startB.endB[,1] > startA.endA.startB.endB[,3])[keep1]
	#	coordinates      <- paste0(ContigAccession,":",StartPosition,"-",EndPosition)
	#	coordinates[take.rv.Comp] <- paste0(ContigAccession,":c",EndPosition,"-",StartPosition)[take.rv.Comp]
	#	targets.mat.temp <- cbind(coordinates,ContigAccession,StartPosition,EndPosition,TargetLength)
	#	result[[i]]      <- targets.mat.temp
	#	setTxtProgressBar(pb,i)
	#}
	#result <- do.call(rbind,result)
	result <- rbind(result.AB,result.BA)
	colnames(result) <- c("coordinates","ContigAccession","StartPosition","EndPosition","TargetLength")
	write.table(result,file=coordinates.output.name,sep="\t",append=F,col.names=T,row.names=F,quote=F)
	result
}
#' SbfI.Seq                              <- REEs::datasets("restriction.enzymes")[datasets("restriction.enzymes")[,"RE.name"]=="SbfI","recognition.sequence"]
#' EcoR1.Seq                             <- REEs::datasets("restriction.enzymes")[datasets("restriction.enzymes")[,"RE.name"]=="EcoRI","recognition.sequence"]
#' Thermophis.baileyi_genome.url         <- REEs::datasets(1)[which(datasets(1)[,"species"]=="Thermophis baileyi"),"genome.url"]
#' proposed.ddRAD.loci.coordinates       <- REEs::proposeLoci.ddRADlike(input.seqs=Thermophis.baileyi_genome.url,output.dir="/ddRAD-like/",recognitionSeqs=c(SbfI.Seq,EcoR1.Seq),lim.lengths=c(900,1000),save.tables=T)





#' In Silico ddRAD version2
#' 
#' Perform in silico ddRAD on a reference genome.
#' This should be closer to actual ddRAD
#' 
#' @param input.seqs A DNAStringSet object or a character string indicating a local filepath or URL path to the input sequences in fasta format.
#' @param output.path Character string indicating where to save the output. If NULL (the default), the output is returned as an object but not saved to a file.
#' @param recognitionSeq A character vector containing the recognition sequence of each restriction enzyme.
#' @param cut.after Numerical vector of length ≥ 1, indicating the position of each cut site within each recognition sequence. A value of zero indicates that the cut site is immediately before the first base of the recognition sequence.
#' @param save.tables Should the tables holding the locations of each recognition sequence on each contig be saved? Default is TRUE but these tables can be big.
#' @param include.antisense Should antisense fragments also be include. Default TRUE. Probably should keep as TRUE for restriction enzymes that produce sticky-ends.
#' @param include.uncut.contigs Should contigs that did not contain any recognition sites be included in the output. Default is TRUE.
#' @param fragment.size.filter Either NULL (default) or a numerical vector of length 2 that defines the minimum and maximum length fragments to keep. If NULL, sequences are not filtered by length.
#' @return A DNAStringSet object containing digested DNA sequences.
#' @export RAD.digestion
RAD.digestion <- function(input.seqs,output.path=NULL,recognitionSeq,cut.after,include.antisense=T,include.uncut.contigs=T,fragment.size.filter=NULL){
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
	names(fa0)           <- gsub(" .+","",names(fa0))
	### Calculate lengths of RE recognition sequences
	RecognitionSeq.length <- nchar(recognitionSeq)
	cut.after.length      <- as.numeric(cut.after)
	### Set up vectors to fill during for loop
	RE.start.pos                <- list(); length(RE.start.pos) <- length(recognitionSeq)
	RE.end.pos                  <- list(); length(RE.end.pos) <- length(recognitionSeq)
	RE.cut.after.pos            <- list(); length(RE.cut.after.pos) <- length(recognitionSeq)
	RE.cut.after.pos.complement <- list(); length(RE.cut.after.pos.complement) <- length(recognitionSeq)
	print("Finding recognition sites")
	pb = txtProgressBar(min = 0, max = length(recognitionSeq), initial = 0)
	for(i in 1:length(recognitionSeq)){
	### Find start positions of each recognition sequence on each contig
		RE.start.pos.temp         <- gregexpr(recognitionSeq[i], fa0)
		names(RE.start.pos.temp)  <- names(fa0)
		RE.start.pos.temp         <- unlist(RE.start.pos.temp)
		### Drop entries with value "-1", which means the recognition site was not found on the contig.
		if(any(RE.start.pos.temp==-1)){
			RE.start.pos.temp <- RE.start.pos.temp[-which(RE.start.pos.temp== -1)]
		}
		RE.start.pos[[i]] <- RE.start.pos.temp
		### Calculate end positions of recognition sequence
		RE.end.pos[[i]]   <- RE.start.pos.temp+RecognitionSeq.length[i]-1
		### Calculate cut-after positions of recognition sequence
		RE.cut.after.pos[[i]] <- (RE.start.pos.temp+cut.after.length[i])-1
		### Calculate cut-after positions of recognition sequence on complement of contigs.
		RE.cut.after.pos.complement[[i]] <- (RE.end.pos[[i]]-cut.after.length[i])
		setTxtProgressBar(pb,i)
	}
	RE.start.pos                <- unlist(RE.start.pos)
	RE.end.pos                  <- unlist(RE.end.pos)
	RE.cut.after.pos            <- unlist(RE.cut.after.pos)
	RE.cut.after.pos.complement <- unlist(RE.cut.after.pos.complement)
	### Create a matrix of contigs, start positions, and end positions for RecognitionSeqA and RecognitionSeqB
	Recognition.data <- cbind(names(RE.start.pos),RE.start.pos,RE.end.pos,RE.cut.after.pos,RE.cut.after.pos.complement)
	### Define column names for Recognition.data
	colnames(Recognition.data) <- c("contig.accession","RE.start.pos","RE.end.pos","RE.cut.after.pos","RE.cut.after.pos.complement")
	### Update the Genbank accession numbers in column 1 of RecognitionSeqA.data and RecognitionSeqB.data.
	Recognition.data[,"contig.accession"] <- gsub("\\.1.+$","\\.1",Recognition.data[,"contig.accession"])
	### Get the list of contigs that have the recognition site
	contigs <- unique(Recognition.data[,"contig.accession"])
	### Lengths of the contigs that have the recognition site
	contigs.with.cutsites.lengths <- width(fa0[contigs])
	### Define vectors to hold contig ranges
	fragments.start       <- list(); length(fragments.start) <- length(contigs)
	fragments.end         <- list(); length(fragments.end)   <- length(contigs)
	accessionList         <- list(); length(accessionList)   <- length(contigs)
	if(include.antisense==T){
		fragments.start.antisense <- list(); length(fragments.start.antisense) <- length(contigs)
		fragments.end.antisense   <- list(); length(fragments.end.antisense)   <- length(contigs)
	}
	print("Calculating fragment ranges")
	pb = txtProgressBar(min = 0, max = length(contigs), initial = 0)
	for(i in 1:length(contigs)){
		contig.temp           <- contigs[i]
		matches.temp          <- which(Recognition.data[,"contig.accession"] == contig.temp)
		RE.cut.after.pos.temp <- sort(RE.cut.after.pos[matches.temp])
		fragments.start[[i]]  <- c(1,RE.cut.after.pos.temp+1)
		fragments.end[[i]]    <- c(RE.cut.after.pos.temp,contigs.with.cutsites.lengths[i])
		accessionList[[i]]    <- rep(contig.temp,(length(matches.temp)+1))
		if(include.antisense==T){
			RE.cut.after.pos.complement.temp <- sort(RE.cut.after.pos.complement[matches.temp])
			fragments.start.antisense[[i]]        <- c(1,(RE.cut.after.pos.complement.temp+1))
			fragments.end.antisense[[i]]          <- c(RE.cut.after.pos.complement.temp,contigs.with.cutsites.lengths[i])
		}
		setTxtProgressBar(pb,i)
	}
	result <- list(); length(result) <- 4
	fragments.start             <- unlist(fragments.start)
	fragments.end               <- unlist(fragments.end)
	accessionList               <- unlist(accessionList)
	subranges                   <- IRanges::IRanges(start=fragments.start ,end=fragments.end,names=accessionList)
	gsubranges                  <- GenomicRanges::GRanges(seqnames=accessionList,ranges=subranges)
	fragments.sense.seqs        <- Biostrings::getSeq(fa0, gsubranges)
	names(fragments.sense.seqs) <- paste0(accessionList,":",fragments.start,"-",fragments.end)
	result[[1]]                 <- fragments.sense.seqs
	if(include.uncut.contigs==T){
		uncut.contigs              <- setdiff(names(fa0),contigs)
		result.sense.uncut.contigs <- fa0[uncut.contigs]
		names(result.sense.uncut.contigs) <- paste0(names(result.sense.uncut.contigs)," uncut")
		result[[2]]                       <- result.sense.uncut.contigs
	}
	if(include.antisense==T){
		fragments.start.antisense       <- unlist(fragments.start.antisense)
		fragments.end.antisense         <- unlist(fragments.end.antisense)
		subranges.antisense             <- IRanges::IRanges(start=fragments.start.antisense ,end=fragments.end.antisense,names=accessionList)
		gsubranges.antisense            <- GenomicRanges::GRanges(seqnames=accessionList,ranges=subranges.antisense)
		fragments.antisense.seqs        <- Biostrings::reverseComplement(Biostrings::getSeq(fa0, gsubranges.antisense))
		names(fragments.antisense.seqs) <- paste0(accessionList,":revComp.of:",fragments.start.antisense,"-",fragments.end.antisense)
		result[[3]]                     <- fragments.antisense.seqs
	}
	if(include.uncut.contigs==T & include.antisense==T){
		result.antisense.uncut.contigs    <- Biostrings::reverseComplement(fa0[uncut.contigs])
		names(result.antisense.uncut.contigs) <- paste0(names(result.antisense.uncut.contigs)," revComp.uncut")
		result[[4]]                           <- result.antisense.uncut.contigs
	}
	result <- collapse.DNAStringSet(result,use.set.names=F)
	if(!is.null(fragment.size.filter)){
		fragment.size.filter <- as.numeric(fragment.size.filter)
		keep.fragments <- which(width(result) > fragment.size.filter[1] & width(result) < fragment.size.filter[2])
		result <- result[keep.fragments]
	}
	if(!is.null(output.path)){
		writeXStringSet(result,output.path)
	}
	result
}
#' @examples
#' SbfI.Seq                              <- datasets("restriction.enzymes")[datasets("restriction.enzymes")[,"RE.name"]=="SbfI","recognition.sequence"]
#' SbfI.cut.after                        <- as.numeric(datasets("restriction.enzymes")[datasets("restriction.enzymes")[,"RE.name"]=="SbfI","cut.after.length"])
#' EcoR1.Seq                             <- datasets("restriction.enzymes")[datasets("restriction.enzymes")[,"RE.name"]=="EcoRI","recognition.sequence"]
#' EcoR1.cut.after                       <- as.numeric(datasets("restriction.enzymes")[datasets("restriction.enzymes")[,"RE.name"]=="EcoRI","cut.after.length"])
#' Thermophis.baileyi_genome.url         <- REEs::datasets(1)[which(datasets(1)[,"species"]=="Thermophis baileyi"),"genome.url"]
#' file.out                              <- "/Users/alyssaleinweber/Documents/SequenceCapture-GitHub/ddRAD-like/Thermophis.baileyi_Digested-genome_CCTGCAGG-GAATTC_900to1000.fas"
#' Thermophis.baileyi_digested.900to1000 <- RAD.digestion(input.seqs=Thermophis.baileyi_genome.url,output.path=file.out,recognitionSeq=c(SbfI.Seq,EcoR1.Seq),cut.after=c(SbfI.cut.after,EcoR1.cut.after),include.antisense=T,include.uncut.contigs=T,fragment.size.filter=c(900,1000))
#' 






