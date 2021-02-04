#' Concatenate DNAStringSets
#' 
#' Concatenate multiple DNA sequence alignments stored in DNAStringSet objects.
#' 
#' @param ... Two or more DNAStringSet objects separated by commas. All taxa must be present in all input DNAStringSets. All sequences in each DNAStringSet must have the same length and are assumed to represent homologous sites.
#' @return A DNAStringSet containing the concatenated sequences in the order in which input DNAStringSets were supplied.
#' @export 
concatenateDNAStringSets <- function(...){
	x              <- list(...)
	x.names        <- lapply(X=x,FUN=names)        ### List of character vectors, each containing the names of the sequences in the corresponding DNAStringSet
	x.names2       <- lapply(X=x.names,FUN=sort)   ### Each character vector of names is sorted alphabetically
	names.matrix   <- do.call(rbind,x.names2)      ### Character matrix. each row is a vector in x.names2
	unique.names   <- unique(names.matrix)         ### 
	if(nrow(unique.names)>1){
		result = print("one or more sequence names not in all input alignments")
	} else {
		x.characters   <- lapply(X=x,FUN=as.character)
		x.mats         <- do.call(cbind,x.characters)
		concat.list    <- apply(X=x.mats,MARGIN=1,FUN=function(y){paste(y,collapse="")})
		result         <- Biostrings::DNAStringSet(concat.list)
	}
	result
}


#' Simple Concatenate DNA Alignment
#' 
#' Used to concatenate multiple DNAStringSet or DNAMultipleAlignment files. Concatenation order is the same as the input argument order.
#' Optionally creates a partition table (option not yet implemented).
#' Compare with function concatenateDNAStringSets.
#' 
#' @param  ... Two or more DNAStringSet objects separated by commas. All taxa must be present in all input DNAStringSets. All sequences in each DNAStringSet must have the same length and are assumed to represent homologous sites.
#' @param  make.partitions.table Logical indicating if a partition table should also be created. Default is FALSE. This argument is not yet implemented.
#' @return A concatenated DNA alignment held in DNAStringSet object.
#' @export 
simple.concatenate.DNA <- function(...,make.partitions.table=F){
	list.of.alignments <- list(...)
	for(i in seq(length(list.of.alignments))){
		names(list.of.alignments)[i]  <-  as.character(sys.call()[[i+1]])
	}
	list.of.loci       <- names(list.of.alignments)

	for(i in 1:length(list.of.alignments)){
		alignment.temp           <- data.table::data.table(as.matrix(list.of.alignments[[i]]),keep.rownames = TRUE)
		locus.name.temp          <- list.of.loci[i]
		col.temp                 <- c(1:(ncol(alignment.temp)-1))
		colnames(alignment.temp) <- c("rn",paste(locus.name.temp,col.temp,sep="."))                          ### redefining column names of data.table so that duplicate colnames wont exist when concatenating

		if(i==1){
			alignment.current   <- alignment.temp
			start.current       <- 1
			end.current         <- (ncol(alignment.current)-1)
			app                 <- FALSE
			rn                  <- alignment.current$rn
			alignment.current   <- apply(X=alignment.current[ ,!"rn"], MARGIN=1, FUN=paste, collapse="")
			alignment.current   <- data.table::data.table(cbind(rn,alignment.current))
		
		} else {
			start.current      <- (end.current+1)
			width.current      <- max(nchar(alignment.current$alignment.current))
			alignment.current1 <- merge(alignment.current, alignment.temp, by="rn", all=TRUE)  ### concatenated data frame
			alignment.current1$alignment.current <- REEs::na.replace(alignment.current1$alignment.current,paste(rep("-",width.current),collapse=""))
			alignment.current2 <- REEs::na.replace(alignment.current1,"-")
			alignment.current  <- alignment.current2
			end.current        <- (start.current+(ncol(alignment.temp2)-2))
			rn                 <- alignment.current$rn
			alignment.current  <- apply(X=alignment.current[ ,!"rn"], MARGIN=1, FUN=paste, collapse="")
			alignment.current  <- data.table::data.table(cbind(rn,alignment.current))
			app                <- TRUE
		}
	}
	alignment.write        <- Biostrings::DNAStringSet(apply(alignment.current[ ,!"rn"], 1, paste, collapse=""))
	names(alignment.write) <- alignment.current$rn
	alignment.write
}