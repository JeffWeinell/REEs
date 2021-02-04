#' Download NCBI Sequences
#' 
#' Download specific regions of a set of sequences from NCBI.
#' 
#' @param outfile Filepath where sequences should be saved.
#' @param input.seqs Either NULL (default but not preferred usage), or a character string indicating either a local filepath or URL path to the input sequences in fasta format. The URL character string is often the preferred method.
#' @param accessionList Character vector of accession numbers.
#' @param startList Numeric vector of starting positions (default = 1)
#' @param endList Numeric vector of end positions.
#' @param strandList Single character or character vector indicating the strand of the sequences (default = "1"). If a single character supplied, then the same strand will be downloaded for all sequences. If length(strandList)>1, then length(strandList) should equal length(accessionList). Set strandList = "1" to download forward direction sequences; strandList = "2" to download reverse complement of sequences.
#' @param db NCBI database to download sequences from (default is "nuccore"). I have not tested this function for other databases. 
#' @param rettype Format to save sequences (default is "fasta"). Alternatively, "gb" (= download sequences in GenBank flatfile format), "native" (= XML format), "acc" (only download accession numbers), "seqid" (download SeqID string), ft" (download feature tables).
#' @param retmode String indicating the mode to download sequences (default = "text"). Don't change from default.
#' @param if.outside.range Ignored if input.seqs is NULL. Otherwise, one of the following character strings: "partial" (the default), "drop", or "error". If partial, then the overlapping portion of sequences are returned if the requested ranges  only partially overlap with the contig. If "error", an error is returned indicating which requested sequences are outside the range of the contig. If "drop", then only requested sequences that are entirely contained within a contig are returned.
#' @param trim.ambiguous Whether or not to trim N strings from the ends of sequences after downloading (default = FALSE). This argument is ignored if input.seqs is NULL.
#' @return A file containing the downloaded sequences.
#' @export get_ncbi_sequences
get_ncbi_sequences <- function(outfile=NULL,input.seqs=NULL,accessionList,startList=1,endList,strandList="1",db="nuccore",rettype="fasta",retmode="text",if.outside.range="partial",trim.ambiguous = FALSE){
	## If input.seqs is NULL (rather than a filepath or URL path), build URLs to fasta formatted sequence for each contig range, and then download each sequence.
	## This is NOT the preferred method.
	if(is.null(input.seqs)){
		URLs <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
			"db=nuccore",
			"&id=",accessionList,
			"&strand=",strandList,
			"&seq_start=",startList,
			"&seq_stop=",endList,
			"&rettype=",rettype,
			"&retmode=",retmode,
			sep=""
		)
		if(is.null(outfile)){
			outfile <- tempfile()
		}
		for(i in 1:length(URLs)){
			temp.URL <- RCurl::getURL(URLs[i])
			write(temp.URL,file=outfile,append=T)
		}
		seqs.scaff.new <- readDNAStringSet(outfile)
	} else {
		## Check if input.seqs is a DNAStringSet object. If true, get and then write subsequences.
		if("DNAStringSet" %in% class(input.seqs)){
			fa0                <- input.seqs
			headers            <- gsub(" .*","",names(fa0))
			scaff.matches.all  <- match(accessionList, headers)
			names(fa0)         <- headers
			fa                 <- fa0[scaff.matches.all]
			delete.subject     <- F
		} else {
			## Check if input.seqs is a local filepath. If true, read input.seqs into R, get and then write subsequences, and then delete the tempfile.
			if(file.exists(input.seqs)){
				### UPDATE THIS ###
				delete.subject       <- F
			} else {
				## Check if input.seqs is a URL path to sequences in fasta format. If true, download seqs to tempfile, read tempfile into R, get and then write subsequences, and then delete the tempfile
				subject.path <- tempfile()
				## sets time limit for downloading files to 1000 seconds
				options(timeout=1000)
				## downloads input.seqs to tempfile
				utils::download.file(url=input.seqs, destfile=subject.path)
				## Get header info for the downloaded sequences
#				fai                  <- Biostrings::fasta.index(subject.path)
#				headers              <- fai$desc
				
				## Read input sequences into R
				fa0                  <- Biostrings::readDNAStringSet(subject.path)
				## Delete all except the accession number from each header string
				headers              <- gsub(" .+","",names(fa0))
				## Which sequences in input.seqs are the ones that we want.
				scaff.matches.all    <- match(accessionList, headers)
				## Read in scaffolds corresponding to accessionList. This takes a minute or two.
				## It may be faster to read in unique values of scaff.matches.all.... still need to figure this out.
#				fa                   <- Biostrings::readDNAStringSet(fai[scaff.matches.all,])
				## Set scaffold sequence names to headers
				names(fa0)           <- headers
				fa                   <- fa0[scaff.matches.all]
				#fa                  <- Biostrings::readDNAStringSet(filepath=subject.path)
				#names(fa)           <- contig.names
				delete.subject       <- T
			}
		}
		start.all          <- apply(X=cbind(startList,endList),MARGIN=1,FUN=min)
		end.all            <- apply(X=cbind(startList,endList),MARGIN=1,FUN=max)
		accession.all      <- accessionList
		if(any(width(fa)   < end.all)){
			### Use this if you want to return an error if any of the requested subsequence ranges are outside of the contig range.
			if(if.outside.range=="error"){
				problem.seqs   <- which(width(fa) < end.all)
				error.message  <- paste0("Requested range(s) of sequence(s) i in ",paste(problem.seqs,collapse=",")," do not overlap (or only partially overlap) contig sequences.")
				stop(error.message)
			}
			### Use this if you want to drop query sequencess if their range is partially outside of the contig range.
			if(if.outside.range=="omit"){
				drop.seqs      <- which(width(fa) < end.all)
				fa             <- fa[-drop.seqs]
				accessionList  <- accessionList[-drop.seqs]
				start.all      <- start.all[-drop.seqs]
				end.all        <- end.all[-drop.seqs]
			}
			if(if.outside.range=="partial"){
			### Use this if you want to include partial sequences if queried sequence range is outside of the contig range.
				partial.seqs          <- which(width(fa) < end.all)
				### Update end.all values if contig length is less than respective end.all value
				end.all[partial.seqs] <- width(fa)[partial.seqs]
			}
		}
		subranges          <- IRanges::IRanges(start=start.all ,end=end.all,names=accessionList)
		gsubranges         <- GenomicRanges::GRanges(seqnames=accessionList,ranges=subranges)
		#### Now need to check if any end.all values are larger than the contig width. Drop sequences for which this is true. ###
#		fa                 <- fa[match(unique(names(fa)),names(fa))]
		seqs.scaff         <- Biostrings::getSeq(fa, gsubranges)
		query.subject.id   <- paste0(accessionList,":",start.all,"-",end.all)
		names(seqs.scaff)  <- query.subject.id
		#close(fa)
		if(trim.ambiguous == T){
			Ns.left.range.all  <- stringr::str_locate(seqs.scaff,"^N+")
			Ns.right.range.all <- stringr::str_locate(seqs.scaff,"N+$")
			### Create a four column matrix showing the start and end positions of left N string (first and second columnd) and right N string (third and fourth columns)
			Ns.ranges.mat   <- cbind(Ns.left.range.all,Ns.right.range.all)
			### Two column matrix holding range of sequences excluding left and right Ns; NA values in first column should be replaced with 1; NA values in second column should be replaced with length of the sequence before any trimming.
			keep.ranges.mat.temp <- cbind((Ns.ranges.mat[,2]+1),(Ns.ranges.mat[,3]-1))
			### Column 1 of keep.ranges.mat.temp, with NA values replaced with 1
			start.keep      <- tidyr::replace_na(keep.ranges.mat.temp[,1],1)
			### Column 2 of keep.ranges.mat.temp, with NA values replaced with respective input sequence length
			end.keep        <- apply(X=cbind(keep.ranges.mat.temp[,2],width(seqs.scaff)),MARGIN=1,FUN=function(row.i){c(NA,row.i)[!is.na(c(NA,row.i))][1]})
			### Two column matrix holding subsequence ranges to keep of seqs.scaff
			#keep.ranges.mat <- cbind(start.keep,end.keep)
			### Amount to shift (right) the start position of target sequence on contig
			start.adjusted <- start.keep-1
			### Amount to shift (left) the end position of target sequence on contig
			end.adjusted   <- width(seqs.scaff)-end.keep
			### Updated values of start.all and end.all after adjusting by start.adjusted and end.adjusted.
			start.all.new  <- start.all+start.adjusted
			end.all.new    <- end.all-end.adjusted
			### Extract subsequence from start.keep to end.keep from respective sequence in seqs.scaff
			subranges.keep  <- IRanges::IRanges(start=start.keep ,end=end.keep,names=query.subject.id)
			gsubranges.keep <- GenomicRanges::GRanges(seqnames=query.subject.id,ranges=subranges.keep)
			seqs.scaff.new  <- Biostrings::getSeq(seqs.scaff, gsubranges.keep)
		} else {
			start.all.new  <- start.all
			end.all.new    <- end.all
			seqs.scaff.new <- seqs.scaff
		}
		
		query.subject.id.new  <- paste0(accessionList,":",start.all.new,"-",end.all.new)
		names(seqs.scaff.new) <- query.subject.id.new
		if(!is.null(outfile)){
			Biostrings::writeXStringSet(x = seqs.scaff.new, filepath=outfile, append=F, format="fasta")
		}
		if(delete.subject){
			file.remove(subject.path)
		}
		#seqs.scaff.new
	}
	result <- seqs.scaff.new
	result
}
