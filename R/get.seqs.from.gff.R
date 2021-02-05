#' Get Sequences from GFF Table
#' 
#' Identical to the function get.exome.from.gff, but used more generally.
#' Note to self: I used this function to extract the MHC loci from the Thamnophis sirtalis genome.
#' 
#' @param input.seqs Character string indicating the local or URL path to input sequences.
#' @param input.gff Either a data table object, or character string  or local path to the filtered gff containing the the loci of interest.
#' @param output.path Path where exomes should be saved (default is NULL).
#' @param additional.ID Filename of a table that cross-references contig names to names in GFF file (default is NULL).
#' @return DNAStringSet with sequences corresponding to unique ranges in input.gff; sequences are written to output.file in fasta format if output.path is not NULL.
#' @export get.seqs.from.gff
get.seqs.from.gff <- function(input.seqs,input.gff,output.path = NULL,additional.ID=NULL) {
	if("character" %in% class(input.gff)){
		filtered.gff1B   <- data.table::fread(input.gff)
	}
	if("data.frame" %in% class(input.gff)){
		filtered.gff1B   <- input.gff
	}
	
	refseq.names     <- unlist(filtered.gff1B[,1])
	if(!is.null(additional.ID)){
		### This is needed because the names in the gff file are not exactly the same as the ones in the genome file; they both link to a common Genbank Accession # though
		ScaffoldKey      <- data.table::fread(additional.ID)
		### Links the names in the gff file to the names of the CDS matches
		alt.scaff.names  <- ScaffoldKey$ScaffoldName[match(refseq.names,ScaffoldKey$RefSeq.ScaffoldAccession)]
		subject.id       <- alt.scaff.names
		subject.start    <- as.numeric(filtered.gff1B$start)
		subject.end      <- as.numeric(filtered.gff1B$end)
		new.names        <- paste0(subject.id,"_",refseq.names,":",subject.start,"-",subject.end)	
	} else {
		subject.id       <- refseq.names
		subject.start    <- as.numeric(filtered.gff1B$start)
		subject.end      <- as.numeric(filtered.gff1B$end)
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
		gr <- as(GenomeInfoDb::seqinfo(fa), "GRanges") ### Check that the seqinfo function isnt actually the one from biofiles package.
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

