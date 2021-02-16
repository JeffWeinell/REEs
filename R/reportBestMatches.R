#' Get Sequences From BLAST Table
#' 
#' Filter a blast hit table of matches to include only the best matches (highest bitscore) to each query.
#' Need to add ability to:
# ' (1) Drop matches if bitscore of top hit is not much better than other matches.
#'  (2) Drop query sequences if the best match has a bitscore below some threshold.
#' 
#' @param input.table Character string with path BLAST hit table, or, a BLAST hit table held as an object of class data.table, data.frame, or matrix.
#' @param output.table.path Character to string with path where output table should be saved; default is NULL.
#' @param remove.subseq.matches Logical with whether or not to drop matches for which the subject sequence is a subsequence of another match. Default is FALSE.
#' @param min.bitscore Number indicating the minimum bitscore required to keep the best match. Matches with a bitscore below this value are not included in the output table. This is useful for removing paralog matches. Default = 50.
#' The reason for using a default min.bitscore of 50 is because values at or greater than this usually indicate that the sequences are homologous (doi: 10.1002/0471250953.bi0301s42).
#' @param min.bitscore.difference Number indicating the minimum difference required between bitscores of the best and second best matches; must be greater than this value to keep matches for the query sequence. Default = 0. This is useful for removing loci with putatitive recent duplicates in the genome.
#' @param blast.method At present, this argument is ignored unless set to "tblastn", in which case the translation frame suffix is removed from query sequence IDs (if present). In the future, I will either delete this argument or this argument will take as input one of the following: NULL (default), or one of "blastn", "blastp", "blastx", "tblastn", "tblastx".
#' @param table.format Number indicating which NCBI table format is applicable to input.table; At present, only format 6 can be used. If more than 12 columns are present only the first 12 are used and assumed to correspond to the columns of format 6.
#' @return A data.table object containing the best matche to each input query in the input blast table. If output.table.path argument is a path (character string), then the function also writes the output as a tab-separated file.
#' @export reportBestMatches
reportBestMatches <- function(input.table, output.table.path=NULL, remove.subseq.matches=F, min.bitscore=50, min.bitscore.difference=0,blast.method=NULL,table.format=6){
	if(table.format!=6){
		stop("At present, only NCBI format 6 tables can be proccessed")
	}
	## Coerce input.table to a data frame object if it is a data.table or matrix object
	if(is(input.table,"data.table") | is(input.table,"matrix")){
		all.matches <- as.data.frame(input.table)
	}
	## If input.table is a character string with the filepath, read the file in as a data.table and then coerce it to data frame.
	if(is(input.table,"character")){
		all.matches       <- as.data.frame(data.table::fread(input=input.table,sep="\t"))
	}
	# Set column names
	colnames(all.matches) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
	##### Set column modes
	# Set which columns should be mode numeric
	numeric.columns <- c(3:12)
	# Set mode to numeric for those columns that should be numeric
	all.matches[, numeric.columns] <- sapply(all.matches[, numeric.columns], as.numeric)
	# Set which columns should be mode character
	character.columns <- c(1:2)
	# Set mode to "character" for the columns indexed in the character.columns vector
	all.matches[, character.columns] <- sapply(all.matches[, character.columns], as.character)
	### In the future, uncomment the next if statement to process hit tables generated from tblastn.
	if(blast.method=="tblastn"){
		all.matches[,"qseqid"] <- mgsub(c("_F1$","_F2$","_F3$","_RC1$","_RC2$","_RC3$"),c("","","","","",""), all.matches[,"qseqid"])
	}
	### Filter matches with bitscore less than min.bitscore
	if(any(all.matches$bitscore < min.bitscore)){
		filtered.matches <- all.matches[-which(all.matches$bitscore < min.bitscore),]
	} else {
		filtered.matches <- all.matches
	}
	# Sort rows, first by decreasing qseqid and then by decreasing bitscore
	matches.ordered    <- filtered.matches[with(filtered.matches, order(filtered.matches[,"qseqid"], filtered.matches[,"bitscore"], decreasing=T)),]
	# Remove rows (matches) if the match is a subsequence of another match.
	if(remove.subseq.matches){
		### argument "split" reguires that qseqid and sseqid columns are mode factor, not mode character.
		# Set mode of qseqid and sseqid columns to factor
		# matches.ordered[, c("qseqid","sseqid")] <- sapply(matches.ordered[, c("qseqid","sseqid")], as.factor)
		# The next line will take a few minutes to complete. It creates a list of data frames. Each data frame includes the rows of matches.ordered for a particular qseqid-sseqid pair.
#		matches.list        <- split(matches.ordered,f=list(matches.ordered$qseqid, matches.ordered$sseqid))
		### Create two column data.frame containing unique pairs from matches.ordered$qseqid and matches.ordered$sseqid
		query.subject.pairs <- dplyr::distinct(matches.ordered[,c("qseqid","sseqid")])
		# The next line will take a few minutes to complete. It creates a list of data frames. Each data frame includes the rows of matches.ordered for a particular qseqid-sseqid pair.
#		matches.list        <- lapply(X=1:nrow(query.subject.pairs),FUN=function(i){ matches.ordered[which(matches.ordered[,"qseqid"] == query.subject.pairs[i,"qseqid"] & matches.ordered[,"sseqid"] == query.subject.pairs[i,"sseqid"]),]})
		matches.list <- list(); length(matches.list) <- nrow(query.subject.pairs)
		for(i in 1:nrow(query.subject.pairs)){
			matches.list[[i]] <- matches.ordered[which(matches.ordered[,"qseqid"] == query.subject.pairs[i,"qseqid"] & matches.ordered[,"sseqid"] == query.subject.pairs[i,"sseqid"]),]
		}
		# May be a good idea to do the next two lines on matches.ordered, rather than on matches.list
		min.list            <- lapply(matches.list,FUN=function(x){apply(X=x[,c("sstart","send")],MARGIN=1,FUN=min)})
		max.list            <- lapply(matches.list,FUN=function(x){apply(X=x[,c("sstart","send")],MARGIN=1,FUN=max)})
		### Create an IRangesList object that is essentially a list of IRanges objects.
		subject.ranges      <- IRanges::IRangesList(start=min.list,end=max.list)
		### Create a HitList object for subject.ranges. The query column of each entry contains the index of ranges that are wholly contained within some other range.
		drop.hits           <- IRanges::findOverlaps(subject.ranges,drop.self=T,type="within")
		### A list of numerical vectors. Each vector contains the rows to drop for the corresponding entry of matches.list
		remove.rows.list    <- lapply(drop.hits,FUN=function(x){unique(as.matrix(x)[,1])})
		### function1 drops the rows in B from a table A as long as B is not empty.
		function1           <- function(A,B){if(length(B)>0){A[-B,]} else{A}}
		### Next line finally performs the filter step
		filtered.matches.ordered.temp  <- t(as.list(mapply(FUN=function1,A=matches.list,B=remove.rows.list)))
		### Next set of lines (until colnames) is getting organizing the data into the correct format
		filtered.matches.ordered.temp2 <- list(); length(filtered.matches.ordered.temp2) <- ncol(filtered.matches.ordered.temp)
		for(i in 1:ncol(filtered.matches.ordered.temp)){
			filtered.matches.ordered.temp2[[i]] <- unlist(filtered.matches.ordered.temp[,i])
		}
		#filtered.matches.ordered           <- data.table::as.data.table(do.call(cbind,filtered.matches.ordered.temp2))
		filtered.matches.ordered           <- as.data.frame(do.call(cbind,filtered.matches.ordered.temp2))
		colnames(filtered.matches.ordered) <- colnames(matches.ordered)
		##### Set column modes
		# Set which columns should be mode numeric
		numeric.columns <- c(3:12)
		# Set mode to numeric for those columns that should be numeric
		filtered.matches.ordered[, numeric.columns] <- sapply(filtered.matches.ordered[, numeric.columns], as.numeric)
		# Set which columns should be mode character
		character.columns <- c(1:2)
		# Set mode to "character" for the columns indexed in the character.columns vector
		filtered.matches.ordered[, character.columns] <- sapply(filtered.matches.ordered[, character.columns], as.character)
	} else {
		filtered.matches.ordered <- matches.ordered
	}
	### Report the best and second best matches.
#	best.matches          <- as.numeric(match(unique(filtered.matches.ordered$qseqid), filtered.matches.ordered$qseqid))
	best.matches          <- as.numeric(match(unique(filtered.matches.ordered$qseqid), filtered.matches.ordered$qseqid))
	best.data             <- filtered.matches.ordered[best.matches,]
	if(min.bitscore.difference!=0){
		without.best.data     <- filtered.matches.ordered[-best.matches,]
		second.best.matches   <- as.numeric(match(unique(filtered.matches.ordered$qseqid), without.best.data$qseqid))
		second.best.data      <- without.best.data[second.best.matches]
		#### If a query does not have a second best match, then a bitscore of zero is used for an dummy second match. There is no consequence for having a single strong match.
		if(any(is.na(second.best.matches))){
			second.best.data$bitscore[which(is.na(second.best.matches))] <- 0
		}
		bitscore.difference             <- as.numeric(best.data$bitscore)-as.numeric(second.best.data$bitscore)
		if(any(bitscore.difference < min.bitscore.difference)){
			best.data <- best.data[-which(bitscore.difference < min.bitscore.difference)]
		}
	}
	if(!is.null(output.table.path)){
		write.table(best.data,file=output.table.path,sep="\t",row.names=F,col.names=T,quote=F)
	}
	# Coerce the result (the best hits table) to a data.table object because it is safe to print.
	result <- data.table::as.data.table(best.data)
	result
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
