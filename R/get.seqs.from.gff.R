#' Get Sequences from GFF Table
#' 
#' Identical to the function get.exome.from.gff, but used more generally.
#' Note to self: I used this function to extract the MHC loci from the Thamnophis sirtalis genome.
#' 
#' @param input.seqs Character string with either the local path or URL to fasta formatted sequences. File can be .gz compressed.
#' @param input.gff One of the following: An object of class data table, data.frame, or character matrix, with columns identical to those used in Generic Feature Format (GFF) tables; or, a character string with local file path to either a GFF file or table with columns identical to those used in GFF files; or, a character string with URL path to a GFF file.
#' @param output.path Character string with path where sequences should be saved (default is NULL).
#' @param additional.ID One of the following: NULL (the default, and usually what you want to use); or, a character string with local filepath to a two-column table with contig names in GFF and sequence file, respectively. Usually, contig names are already the same in GFF and sequence files, but this argument can be useful if you renamed contigs in either the GFF or sequence file.
#' @return DNAStringSet object with sequences corresponding to unique ranges in input.gff; sequences are written to output.file in fasta format if output.path is a character string.
#' @export get.seqs.from.gff
get.seqs.from.gff <- function(input.seqs,input.gff,output.path = NULL,additional.ID=NULL) {
	#if("character" %in% class(input.gff)){
	#	### Read the GFF as a 
	#	filtered.gff1B   <- data.table::fread(input.gff)
	#	### Now change to 
	#}
	#if("data.frame" %in% class(input.gff)){
	#	filtered.gff1B   <- input.gff
	#}
	# Set x to be input.gff because load.gff has an argument also called input.gff
	x <- input.gff
	# Load the GFF as a data.table object and then coerce to a data.frame
	filtered.gff1B <- as.data.frame(load.gff(x))
	##### Set column modes
	# Identify which columns are named "start" or "end", because there are the columnd that should be mode "numeric"
	numeric.columns <- which(colnames(filtered.gff1B) %in% c("start","end"))
	# Set mode to numeric for those columns that should be numeric
	filtered.gff1B[, numeric.columns] <- sapply(filtered.gff1B[, numeric.columns], as.numeric)
	# Get column indices for all of the columns other than the columns with names "start" and "end".
	character.columns <- which(!(colnames(filtered.gff1B) %in% c("start","end")))
	# Set mode to "character" for the columns indexed in the character.columns vector
	filtered.gff1B[, character.columns] <- sapply(filtered.gff1B[, character.columns], as.character)
	# Get contig names. These should be in the first column, which is usually named "seqid"
	refseq.names   <- unlist(filtered.gff1B[,1])
	### Load additional.ID matrix if contig names are not the same in input.seqs and input.gff
	# I havent check this in a while so this may not work as expected.
	if(!is.null(additional.ID)){
		### This is needed because the names in the gff file are not exactly the same as the ones in the genome file; they both link to a common Genbank Accession # though
		ScaffoldKey      <- data.table::fread(additional.ID)
		### Links the names in the gff file to the names of the CDS matches
		alt.scaff.names  <- ScaffoldKey$ScaffoldName[match(refseq.names,ScaffoldKey$RefSeq.ScaffoldAccession)]
		subject.id       <- alt.scaff.names
		subject.start    <- filtered.gff1B$start
		subject.end      <- filtered.gff1B$end
		new.names        <- paste0(subject.id,"_",refseq.names,":",subject.start,"-",subject.end)
	} else {
		subject.id       <- refseq.names
		subject.start    <- filtered.gff1B$start
		subject.end      <- filtered.gff1B$end
		new.names        <- paste0(refseq.names,":",subject.start,"-",subject.end)
	}
	if("DNAStringSet" %in% class(input.seqs)){
		subject.path <- tempfile()
		# fa               <- input.seqs
		# headers          <- names(input.seqs)
		# contig.names     <- mat.strsplit(headers)[,1]
		# names(fa)        <- contig.names
		Biostrings::writeXStringSet(x = subject, filepath=subject.path, append=F, format="fasta")
		delete.subject <- T
	} else {
		if(file.exists(input.seqs)){
			subject.path   <- input.seqs
			delete.subject <- F
		} else {
			subject.path <- tempfile()
			# Sets time limit for downloading files to 1000 seconds
			options(timeout=1000)
			conn <- utils::download.file(url=input.seqs, destfile=subject.path)
			delete.subject <- T
		}
	}
	input.seqs.path <- subject.path
	if(summary(file(input.seqs.path))$class != "gzfile"){
		### Create an index of file 'foo.fasta'; this avoids having to actually copy or move the file to a new directory
		Rsamtools::indexFa(input.seqs)
		fa <- Rsamtools::FaFile(input.seqs)
		gr <- as(GenomeInfoDb::seqinfo(fa), "GRanges")
		contig.names <- names(gr)
	} else {
		### Extract fasta headers (description lines) and then contig names (IDs) from the input seqs file (usually a genome file)
		headers        <- Biostrings::fasta.index(input.seqs.path)$desc
		### Next few lines update headers so that the names follow the same rules for seq_ids that were required by makeblastdb
		#headers <- gsub(" ","_",headers)
		#contig.names <- substring(headers,first=1,last=50)
		# The next line, which is commented out, would work if input.seqs is always an NCBI genome, and
		# contig.names         <- mat.strsplit(headers,split="_")[,1]
		### It takes about 30 seconds to load the entire genome into R on my laptop. This is the longest step of this function.
		fa             <- Biostrings::readDNAStringSet(filepath=input.seqs.path)
		
		# contig.names  <- mat.strsplit(headers)[,1]
		contig.names   <- gsub(" .+","",headers)
		names(fa)      <- contig.names
	}
	scaff.matches.all  <- match(subject.id, contig.names)
	start.all          <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=min)
	end.all            <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=max)
	subranges          <- IRanges::IRanges(start=start.all ,end=end.all,names=subject.id)
	gsubranges         <- GenomicRanges::GRanges(seqnames=subject.id,ranges=subranges)
	output.seqs        <- Biostrings::getSeq(fa, gsubranges)
	names(output.seqs) <- new.names
	output.seqs        <- unique(output.seqs)
	if(!is.null(output.path)){
		Biostrings::writeXStringSet(x = output.seqs, filepath=output.path, append=F, format="fasta")
	}
	if(delete.subject){
		file.remove(subject.path)
	}
	output.seqs
}
#' @examples
#' ### Load GFF table from NCBI repository.
#' Thamnophis.sirtalis_GFF.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"
#' Thamnophis.sirtalis_GFF     <- load.gff(input=Thamnophis.sirtalis_GFF.url,local=F)
#' 
#' # Filter Thamnophis.sirtalis_GFF to only include CDS features with length at least 120bp
#' Thamnophis.sirtalis_GFF_CDS_longer120bp <- filter.gff(input.gff=Thamnophis.sirtalis_GFF,feature.type="CDS",min.length=120)
#' 
#' ### Use get.seqs.from.gff to extract the sequences for the loci in the filtered GFF.
#' Thamnophis.sirtalis_genome.path <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz"
#' Thamnophis.sirtalis_exome <- get.seqs.from.gff(input.seqs=Thamnophis.sirtalis_genome.path,input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp)

