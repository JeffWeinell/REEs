#' Concatenate All Alignments in Directory
#' 
#' Concatenate all phylip-formatted DNA or AA alignments that are stored in separate input files in the same input directory, and optionally create a partition file.
#' 
#' @param alignment.folder Path to the directory containing the input set of alignments that should be concatenated.
#' @param type Type of sequences (default = "DNA"). Alternatively "AA" for Amino Acid sequences.
#' @param seqs.format Format of input sequences (default = "phylip"). No other formats currently implemented.
#' @param alignment.output Where to save the concatenated alignment (default value NA does not save the concatenated alignment).
#' @param partition.output Where to save the partition file for the concatenated alignment (default value NA does not create a partition file).
#' @param show.progress Print progress to screen (default TRUE)
#' @param istart Which alignment to start at (default is NA, meaning all alignments in alignment.folder should be concatenated)
#' @param iend Which alignment to end at (default is NA, meaning all alignments in alignment.folder should be concatenated)
#' @param only.variable.sites Should invarient sites be removed? Default is FALSE.
#' @return DNAStringSet or AAStringSet object containing the concatenated sequences in the order of list.files(alignment.folder).
#' @export concatenate.alignments
concatenate.alignments    <- function(alignment.folder,type="DNA",seqs.format="phylip",alignment.output=NA,partition.output=NA,show.progress=T,istart=NA,iend=NA,only.variable.sites=F){
	alignment.filenames   <- list.files(alignment.folder,full.names=T)
	#alignment.shortnames  <- gsub("\\..+","",list.files(alignment.folder))
	alignment.shortnames  <- nameFromPath(alignment.filenames,include.extension=F)

	if(is.na(istart)){
		istart <- 1
	}
	if(is.na(iend)){
		iend <- length(alignment.filenames)
	}
	
	for(i in istart:iend){
		locus.name.temp <- alignment.shortnames[i]
		if(show.progress==T & is.wholenumber(i/50)){
			print(paste0(i," of ",length(alignment.filenames)))
		}
		if(type=="DNA"){
			alignment.temp   <- Biostrings::readDNAMultipleAlignment(alignment.filenames[i],format=seqs.format)
		}
		if(type=="AA"){
			alignment.temp   <- Biostrings::readAAMultipleAlignment(alignment.filenames[i],format="fasta")
		}
		if(all(names(Biostrings::unmasked(alignment.temp)) %in% c(1:nrow(alignment.temp)))){
			next
		} else {
			if(only.variable.sites==T){
				if(type=="DNA"){
					alignment.temp2      <- filter.alignment(alignment.temp,min.allele.freqs.dna=c(1,1,0,0))
				}
				if(type=="AA"){
					alignment.temp2      <- filter.alignment(alignment.temp,min.allele.freqs.aa=c(c(1,1),rep(0,19)))
				}
			}
			if(only.variable.sites==F){
				alignment.temp2      <- filter.alignment(alignment.temp)
			}
			
		}
		if(ncol(alignment.temp2)==0){
			next
		} else {
			alignment.temp2  <- data.table::data.table(as.matrix(alignment.temp2),keep.rownames = TRUE)    ### Coerces alignment into data.table
			col.temp  <- c(1:(ncol(alignment.temp2)-1))                                                    ### Subtract 1 because first column of data.table is sample name
			### Next line is used for SnakeCap loci but wont harm other data such as Sanger data...
			#key.temp  <- gsub("WeinellEntry","WE",locus.name.temp)                                        ### Shorter version of SnakeCap locus name but wont be generalizable across datasets
			#colnames(alignment.temp2) <- c("rn",paste(key.temp,col.temp,sep="."))                         ### Redefining column names of data.table so that duplicate colnames wont exist when concatenating
			colnames(alignment.temp2) <- c("rn",paste0("L",i,".",col.temp))                                ### Redefining column names of data.table so that duplicate colnames wont exist when concatenating
		}
		if(i==istart){
			alignment.current   <- alignment.temp2
			start.current       <- 1
			end.current         <- (ncol(alignment.current)-1) 
			app                 <- FALSE
			rn                  <- alignment.current$rn
			alignment.current   <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
			alignment.current   <- data.table::data.table(cbind(rn,alignment.current))
		
		} else {
			start.current      <- (end.current+1)
			width.current      <- max(nchar(alignment.current$alignment.current))
			alignment.current1 <- merge(alignment.current, alignment.temp2, by="rn", all=TRUE)  ### concatenated data frame
			#alignment.current2  <- na.replace(alignment.current1$alignment.current,rep("-",width.current))
			
			alignment.current1$alignment.current <- REEs::na.replace(x=alignment.current1$alignment.current,replace=paste(rep("-",width.current),collapse=""))  #### Indivduals without a locus are included with gaps.
			#alignment.current <- na.replace(merge(alignment.current, alignment.temp2, by="rn", all=TRUE),"-")  #### Concatenated data frame
			
			alignment.current2 <- REEs::na.replace(x=alignment.current1,replace="-")
			alignment.current  <- alignment.current2
			end.current        <- (start.current+(ncol(alignment.temp2)-2))
			rn                 <- alignment.current$rn
			alignment.current  <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
			alignment.current  <- data.table::data.table(cbind(rn,alignment.current))
			app                <- TRUE
		}
		if(type=="DNA"){
			test.rows <- Biostrings::DNAStringSet(apply(alignment.current[ ,!"rn"], 1, paste, collapse=""))
		}
		if(type=="AA"){
			test.rows <- Biostrings::AAStringSet(apply(alignment.current[ ,!"rn"], 1, paste, collapse=""))
		}
		### Sanity check that sequences are all the same length. Throws an error they are not.
		if(!all(width(test.rows)==width(test.rows[1]))){
			print(paste(i, "error",sep=" "))
		}
		### Line to write in the partition file.
		text.current = paste(type,", ",locus.name.temp,"_VarSites = ",start.current,"-",end.current,sep="")
		if(!is.na(partition.output)){
			write(text.current,file=partition.output,append=app)
		}
	}
	alignment.temp3        <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
	if(type=="DNA"){
		alignment.write    <- Biostrings::DNAStringSet(alignment.temp3)
	}
	if(type=="AA"){
		alignment.write    <- Biostrings::AAStringSet(alignment.temp3)
	}
	names(alignment.write) <- alignment.current$rn
	if(!is.na(alignment.output)){
		Biostrings::writeXStringSet(x=alignment.write, filepath=alignment.output, append = F)
	}
	alignment.write
}
