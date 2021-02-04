#' Get Sequences From BLAST Table
#' 
#' Filter a blast hit table of matches to include only the best matches (highest bitscore) to each query.
#' Need to add ability to:
# ' (1) Drop matches if bitscore of top hit is not much better than other matches.
#'  (2) Drop query sequences if the best match has a bitscore below some threshold.
#' 
#' @param input.table Path to the input hit table of matches, or, a data.table object.
#' @param output.table.path Where to save the output table of best matches (default is NULL).
#' @param remove.subseq.matches Whether or not to drop matches for which the subject sequence is a subsequence of another match. Default is TRUE.
#' @param min.bitscore Minimum bitscore required for the best match. Matches with a bitscore below this value are not included in the output table. This is useful for removing paralog matches. Default = 50.
#' The reason for using a default min.bitscore of 50 is because values at or greater than this usually indicate that the sequences are homologous (doi: 10.1002/0471250953.bi0301s42).
#' @param min.bitscore.difference The difference between the best and second best matches' bitscore must be greater than this value to keep matches for the query sequence. Default = 0. This is useful for removing loci with putatitive recent duplicates in the genome.
#' @return A data.table object containing the best matche to each input query in the input blast table. If output.table.path argument is a path (character string), then the function also writes the output as a tab-separated file.
#' @export reportBestMatches
reportBestMatches <- function(input.table,output.table.path=NULL,remove.subseq.matches=T, min.bitscore=50,min.bitscore.difference=0){
	if("data.frame" %in% class(input.table)){
		all.matches <- data.table::as.data.table(input.table)
	}
	if("character" %in% class(input.table)){
		all.matches       <- data.table::fread(input=input.table,sep="\t")
	}
	colnames(all.matches) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
	
	#print("check 0")
	### Filter matches with bitscore less than min.bitscore
	if(any(all.matches$bitscore < min.bitscore)){
		filtered.matches <- all.matches[-which(all.matches$bitscore < min.bitscore)]
	} else {
		filtered.matches <- all.matches
	}
	
#	# New method to sort rows. Need to verify that this produced the same result as the old method. This should be slightly faster though.
#	all.matches.ordered  <- all.matches[with(all.matches, order(qseqid, bitscore, decreasing=T)),]

	# New method to sort rows. Need to verify that this produced the same result as the old method. This should be slightly faster though.
	matches.ordered  <- filtered.matches[with(filtered.matches, order(qseqid, bitscore, decreasing=T)),]
#	query.subject.sstart.send <- cbind(matches.ordered$qseqid, matches.ordered$sseqid, matches.ordered$sstart, matches.ordered$send)
	#print("check 1")
	# Remove rows (matches) if the match is a subsequence of another match.
	if(remove.subseq.matches){
		#sstart.send    <- cbind(matches.ordered$sstart, matches.ordered$send)
		#min.vals       <- apply(X=sstart.send,MARGIN=1,FUN=min)
		#max.vals       <- apply(X=sstart.send,MARGIN=1,FUN=max)
		#subject.ranges <- IRanges(start=min.vals, end= max.vals, names= paste(matches.ordered$qseqid, matches.ordered$sseqid,sep="_"))
		###
		#query.subject       <- paste(matches.ordered$qseqid,matches.ordered$sseqid,sep=".")
		matches.list        <- split(matches.ordered,by=c("qseqid","sseqid"))
		min.list            <- lapply(matches.list,FUN=function(x){apply(X=x[,c("sstart","send")],MARGIN=1,FUN=min)})
		max.list            <- lapply(matches.list,FUN=function(x){apply(X=x[,c("sstart","send")],MARGIN=1,FUN=max)})
		### Create an IRangesList object that is essentially a list of IRanges objects.
		subject.ranges      <- IRanges::IRangesList(start=min.list,end=max.list)
		### Create a HitList object for subject.ranges. The query column of each entry contains the index of ranges that are wholly contained within some other range.
		drop.hits           <- IRanges::findOverlaps(subject.ranges,drop.self=T,type="within")
		### A list of numerical vectors. Each vector contains the rows to drop for the corresponding entry of matches.list
		remove.rows.list    <- lapply(drop.hits,FUN=function(x){unique(as.matrix(x)[,1])})
		### function1 drops the rows in B from a table A as long as B is not empty.
		function1 <- function(A,B){if(length(B)>0){A[-B,]} else{A}}
		### Next line finally performs the filter step
		filtered.matches.ordered.temp  <- t(as.list(mapply(FUN=function1,A=matches.list,B=remove.rows.list)))
		### Next set of lines (until colnames) is getting organizing the data into the correct format
		filtered.matches.ordered.temp2 <- list(); length(filtered.matches.ordered.temp2) <- ncol(filtered.matches.ordered.temp)
		for(i in 1:ncol(filtered.matches.ordered.temp)){
			filtered.matches.ordered.temp2[[i]] <- unlist(filtered.matches.ordered.temp[,i])
		}
		filtered.matches.ordered           <- data.table::as.data.table(do.call(cbind,filtered.matches.ordered.temp2))
		colnames(filtered.matches.ordered) <- colnames(matches.ordered)
		#filtered.matches.ordered      <- data.table::data.table(unlist(filtered.matches.ordered.temp[,1]),unlist(filtered.matches.ordered.temp[,2]),unlist(filtered.matches.ordered.temp[,3]))
		#remove.ID           <- unlist(remove.rows.list)
		#query.names         <- gsub("\\.1[0-9]+$","\\.1",names(remove.ID))
		#query.vals          <- as.numeric(remove.ID)
		#names.vals.mat      <- data.frame(query.names,query.vals)
		
		#rows.indices.list   <- lapply(query.names,FUN=function(x){which(query.subject %in% x)})
		#remove.rows.indices <- mapply(FUN=function(A,B){A[B]},A=rows.indices.list,B=query.vals)
		
		### The next line takes a while to run. It holds the indices of the rows to drop of matches.ordered.
		# remove.rows.indices <- mapply(FUN=function(A,B){grep(A,query.subject)[B]},A=query.names,B=query.vals)

		### Possible alternative to remove.rows.indices
		# remove.rows.indices <- apply(X=names.vals.mat,MARGIN=1,FUN=function(A){grep(A[1],query.subject)[A[2]]})
		#test.ranges  <- IRanges(start=min.list[[1]],end=max.list[[1]])
		#if(length(findOverlaps(test.ranges,drop.self=T,type="within"))>0){
		#	test.ranges2 <- test.ranges[-unique(as.matrix(findOverlaps(test.ranges,drop.self=T,type="within"))[,1])]
		#} else {
		#	test.ranges2 <- test.ranges
		#}

		### Unique combination of the query sequences and subject contigs in the hit table
#		distinct.mat <- dplyr::distinct(filtered.matches[,c(1:2)]) 
#		result       <- list(); length(result)=nrow(matches.ordered)
#		print("filtering subsequences")
#		pb = txtProgressBar(min = 1, max = nrow(distinct.mat), initial = 0)
#		for(i in 1:nrow(distinct.mat)){
#			matches.temp       <- matches.ordered[which(matches.ordered$qseqid==distinct.mat$qseqid[i] & matches.ordered$sseqid== distinct.mat$sseqid[i]),]
#			matches.temp.mins  <- apply(X=matches.temp[,c("sstart","send")],MARGIN=1,FUN=min)
#			matches.temp.maxes <- apply(X=matches.temp[,c("sstart","send")],MARGIN=1,FUN=max)
#			for(k in 1:nrow(matches.temp)){
#				#r1.temp                <- c(matches.temp.mins[k],matches.temp.maxes[k])
#				result.index           <- which(matches.ordered$qseqid==matches.temp$qseqid[k] & matches.ordered$sseqid== matches.temp$sseqid[k] & matches.ordered$sstart == matches.temp$sstart[k] & matches.ordered$send == matches.temp$send[k])
#				#result[result.index]   <- any(apply(X=Ophiophagus.hannah.50hits[-11,c("sstart","send")],MARGIN=1,FUN=is.subset,r1=Ophiophagus.hannah.50hits[11,c("sstart","send")]))
#				result[result.index]   <- any(apply(X=matches.temp[-k,c("sstart","send")],MARGIN=1,FUN=is.subset,r1=matches.temp[k,c("sstart","send")]))
#			}
#		setTxtProgressBar(pb,i)
#		}
#		result <- unlist(result)
#		if(any(result)){
#			matches.ordered <- matches.ordered[!result,]
#		}
		#filtered.matches.ordered <- matches.ordered[-remove.rows.indices,]
	} else {
		filtered.matches.ordered <- matches.ordered
	}
#	#print("check 2")
#	### Filter matches with bitscore less than min.bitscore
#	if(any(all.matches.ordered$bitscore < min.bitscore)){
#		filtered.matches.ordered <- all.matches.ordered[-which(all.matches.ordered$bitscore < min.bitscore)]
#	} else {
#		filtered.matches.ordered <- all.matches.ordered
#	}
	
	#print("check 3")
	### Report the best and second best matches.
	best.matches          <- as.numeric(match(unique(filtered.matches.ordered$qseqid), filtered.matches.ordered$qseqid))
	best.data             <- filtered.matches.ordered[best.matches,]
	if(min.bitscore.difference!=0){
		without.best.data     <- filtered.matches.ordered[-best.matches,]
		second.best.matches   <- as.numeric(match(unique(filtered.matches.ordered$qseqid), without.best.data$qseqid))
		second.best.data      <- without.best.data[second.best.matches]
		#print("check 4")
		#### If a query does not have a second best match, then a bitscore of zero is used for an dummy second match. There is no consequence for having a single strong match.
		if(any(is.na(second.best.matches))){
			second.best.data$bitscore[which(is.na(second.best.matches))] <- 0
		}
		#print("check 5")
		bitscore.difference             <- as.numeric(best.data$bitscore)-as.numeric(second.best.data$bitscore)
		if(any(bitscore.difference < min.bitscore.difference)){
			best.data <- best.data[-which(bitscore.difference < min.bitscore.difference)]
		}
	}
	#print("check 6")
	if(!is.null(output.table.path)){
		write.table(best.data,file=output.table.path,sep="\t",row.names=F,col.names=T,quote=F)
	}
	best.data
}
#' @examples
#' ### Load GFF table from NCBI repository.
#' Thamnophis.sirtalis_GFF.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"
#' Thamnophis.sirtalis_GFF     <- load.gff(input=Thamnophis.sirtalis_GFF.url,local=F)
#' 
#' # Filter Thamnophis.sirtalis_GFF to only include CDS features with length at least 120bp
#' Thamnophis.sirtalis_GFF_CDS_longer120bp <- filter.gff(input.gff=Thamnophis.sirtalis_GFF,feature.type="CDS",min.length=120)
#' 
#' ### Use get.exome.from.annotationTable to extract the sequences for the loci in the filtered GFF
#' Thamnophis.sirtalis_genome.path <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz"
#' Thamnophis.sirtalis_exome <- get.seqs.from.gff(input.seqs=Thamnophis.sirtalis_genome.path,input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp)
#' 
#' ### Use tblastx to return up to 50 matches for each of exon (only first two exons in this example) of Thamnophis exome in Crotalus horridus genome
#' test.query    <- Thamnophis.sirtalis_exome[1:2]
#' test.subject  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/625/485/GCA_001625485.1_ASM162548v1/GCA_001625485.1_ASM162548v1_genomic.fna.gz"
#' test.50hits   <- blast(method="tblastx",subject=test.subject,query=test.query)
#' 
#' ### Return the best match to each query sequence
#' best.hits <- reportBestMatches(test.50hits)
