#' @title Get Sequences From BLAST Table
#' 
#' Extracts the sequences defined in an input table (a blast table with NCBI format 6) from an input genome. The input genome should be the same genome used as the subject in the blast table.
#' The quality/completeness of the output sequence dataset depends on the quality/completeness of the blast search.
#' This function supercedes the functions get.exome.from.blastTable and get.UCEs.from.blastTable.
#' 
#' @param input.blastTable Either a character string with the local filepath to a BLAST table (BLAST output format 6), or such a table held as a data table, data frame, or character matrix object.
#' @param input.seqs Either a character string with the filepath or URL path to the input sequences (which can be .gz compressed), or, a DNAStringSet object holding the sequences.
#' @param output.path Character string with path (including filename) to where extracted sequences should be saved.
#' @return A DNAStringSet object with the extracted sequences. Sequences are also optionally written to output.path.
#' @export get.seqs.from.blastTable
get.seqs.from.blastTable <- function(input.blastTable,input.seqs,output.path=NULL){
	#if("data.frame" %in% class(input.blastTable)){
	#	data.best   <- data.table::as.data.table(input.blastTable)
	#}
	#if("character" %in% class(input.blastTable)){
	#	data.best   <- data.table::fread(file=input.blastTable,sep=",",header=T)
	#}
	## Coerce input.blastTable to a data frame object if it is a data.table or character matrix object
	if(is(input.blastTable,"data.table") | is(input.blastTable,"matrix")){
		data.best <- as.data.frame(input.blastTable)
	}
	## If input.blastTable is a character string with the filepath, read the file as a data.table and then coerce it to data frame.
	if(is(input.blastTable,"character")){
		data.best <- as.data.frame(data.table::fread(input=input.blastTable,sep="\t"))
	}
	## Now set filenames and 
	colnames(data.best) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
	##### Set column modes
	# Set which columns should be mode numeric
	numeric.columns <- c(3:12)
	# Set mode to numeric for those columns that should be numeric
	data.best[, numeric.columns] <- sapply(data.best[, numeric.columns], as.numeric)
	# Set which columns should be mode character
	character.columns <- c(1:2)
	# Set mode to "character" for the columns indexed in the character.columns vector
	data.best[, character.columns] <- sapply(data.best[, character.columns], as.character)
	
	query.subject.id    <- paste0(as.character(data.best$qseqid),"_Subject=",as.character(data.best$sseqid),"_",as.character(data.best$sstart),"_",as.character(data.best$send))
	subject.id          <- gsub(" .+","",as.character(data.best$sseqid))
	subject.start       <- data.best$sstart
	subject.end         <- data.best$send
	if("DNAStringSet" %in% class(input.seqs)){
		headers             <- gsub(" .+","",names(input.seqs))
		names(input.seqs)   <- headers
		scaff.matches.all   <- match(subject.id, headers)
		fa                  <- input.seqs[scaff.matches.all]
		delete.subject      <- F
	} else {
		if(file.exists(input.seqs)){
			if(summary(file(input.seqs))$class != "gzfile"){
				### Create an index of file 'foo.fasta'; this avoids having to actually copy or move the file to a new directory
				Rsamtools::indexFa(input.seqs)
				fa <- Rsamtools::FaFile(input.seqs)
				### Check that the seqinfo function isnt actually the one from biofiles package.
				gr <- as(GenomeInfoDb::seqinfo(fa), "GRanges")
				contig.names <- names(gr)
				headers      <- gsub(" .+","",names(gr))
			} else {
				### Extract fasta headers (description lines) and then contig names (IDs) from the genome file
				fa0                <- Biostrings::readDNAStringSet(input.seqs)
				headers            <- gsub(" .+","",names(fa0))
				names(fa0)         <- headers
				scaff.matches.all  <- match(subject.id, headers)
				fa                 <- fa0[scaff.matches.all]
				# The next line, which is commented out, would work if input.seqs is always an NCBI genome, and
				# contig.names         <- mat.strsplit(headers,split="_")[,1]
				### It takes about 30 seconds to load the entire genome into R on my laptop. This is the longest step.
				#fa         <- Biostrings::readDNAStringSet(filepath=input.seqs)
				#names(fa)  <- contig.names
			}
			delete.subject <- F
		} else {
			subject.path <- tempfile()
			# Sets time limit for downloading files to 1000 seconds
			options(timeout=1000)
			conn <- utils::download.file(url=input.seqs, destfile=subject.path)
			## Read input sequences into R
			fa0                  <- Biostrings::readDNAStringSet(subject.path)
			## Close connection to downloaded file to prevent warning about unused connections.
			# Close(conn)
			## Update headers to only include NCBI accession
			headers              <- gsub(" .+","",names(fa0))
			names(fa0)           <- headers
			#headers             <- fai$desc
			#headers             <- mgsub(c(" ",","),c("_","_"),headers)
			#contig.names        <- substring(headers,first=1,last=50)
			#contig.names        <- gsub(" .+","",headers)
			scaff.matches.all    <- match(subject.id, headers)
			# scaff.matches.all  <- pmatch(subject.id, contig.names)
			# scaff.matches.all  <- grep(subject.id, contig.names)
			#fa                  <- Biostrings::readDNAStringSet(fai[scaff.matches.all,])
			#names(fa)           <- contig.names[scaff.matches.all]
			#names(fa)           <- headers
			#fa                  <- Biostrings::readDNAStringSet(filepath=subject.path)
			fa                   <- fa0[scaff.matches.all]
			#names(fa)           <- contig.names
			delete.subject       <- T
		}
	}
	
	#Rsamtools::indexFa(input.seqs)
	#fa <- Rsamtools::FaFile(input.seqs)
	#gr <- as(GenomeInfoDb::seqinfo(fa), "GRanges")
	#scaff.matches.all  <- match(subject.id, contig.names)
	start.all          <- as.numeric(apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=min))
	end.all            <- as.numeric(apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=max))
	subranges          <- IRanges::IRanges(start=start.all ,end=end.all,names=subject.id)
	gsubranges         <- GenomicRanges::GRanges(seqnames=subject.id,ranges=subranges)
	seqs.scaff         <- Biostrings::getSeq(fa, gsubranges)
	names(seqs.scaff)  <- query.subject.id
	#close(fa)
	if(!is.null(output.path)){
		Biostrings::writeXStringSet(x = seqs.scaff, filepath=output.path, append=F, format="fasta")
	}
	if(delete.subject){
		file.remove(subject.path)
	}
	seqs.scaff
}
#' @examples
#' 
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
#' ### Use tblastx to return up to 50 matches for each of query sequence in a subject sequence database
#' # Define the query sequences. In this example, query sequences are the first sequences in the Thamnophis sirtalis exome DNAStringSet object.
#' test.query    <- Thamnophis.sirtalis_exome[1:2]
#' 
#' # Define which sequences to blast against (subject sequences). In this case, the we provide the URL path to the Crotalus horridus genome (the subject sequences).
#' Crotalus.horridus.genome_url  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/625/485/GCA_001625485.1_ASM162548v1/GCA_001625485.1_ASM162548v1_genomic.fna.gz"
#'
#' # URLs to other example genomes (including the Crotalus.horridus.genome_url just defined) are included in this package and can be accessed with the datasets() function:
#' datasets(1) ### Prints a list of URLs to whole genomes of 11 snake species and 3 lizard species.
#' 
#' Anolis.carolinensis.genome_url           <- datasets(1)[2,2]
#' Gekko.japonicus.genome_url               <- datasets(1)[3,2]
#' Pagona.vitticeps.genome_url              <- datasets(1)[4,2]
#' Crotalus.mitchellii.genome_url           <- datasets(1)[6,2]
#' Ophiophagus.hannah.genome_url            <- datasets(1)[7,2]
#' Pantherophis.guttatus.genome_url         <- datasets(1)[8,2]
#' Protobothrops.mucrosquamatus.genome_url  <- datasets(1)[9,2]
#' Python.bivittatus.genome_url             <- datasets(1)[10,2]
#' Viperus.berus.genome_url                 <- datasets(1)[11,2]
#'
#' # Running the blast function.
#' Anolis.carolinensis.50hits          <- blast(method="tblastx",subject=Anolis.carolinensis.genome_url,query=test.query)
#' Gekko.japonicus.50hits              <- blast(method="tblastx",subject=Gekko.japonicus.genome_url,query=test.query)
#' Pagona.vitticeps.50hits             <- blast(method="tblastx",subject=Pagona.vitticeps.genome_url,query=test.query)
#' Crotalus.horridus.50hits            <- blast(method="tblastx",subject=Crotalus.horridus.genome_url,query=test.query)
#' Crotalus.mitchellii.50hits          <- blast(method="tblastx",subject=Crotalus.mitchellii.genome_url,query=test.query)
#' Ophiophagus.hannah.50hits           <- blast(method="tblastx",subject=Ophiophagus.hannah.genome_url,query=test.query)
#' Pantherophis.guttatus.50hits        <- blast(method="tblastx",subject=Pantherophis.guttatus.genome_url,query=test.query)
#' Protobothrops.mucrosquamatus.50hits <- blast(method="tblastx",subject=Protobothrops.mucrosquamatus.genome_url,query=test.query)
#' Python.bivittatus.50hits            <- blast(method="tblastx",subject=Python.bivittatus.genome_url,query=test.query)
#' Viperus.berus.50hits                <- blast(method="tblastx",subject=Viperus.berus.genome_url,query=test.query)

#' ### Filter results to include only the best match for each query sequence
#' best.hits.Crotalus.horridus      <- reportBestMatches(Crotalus.horridus.50hits)
#' best.hits.Crotalus.mitchellii    <- reportBestMatches(Crotalus.mitchellii.50hits)
#' best.hits.Ophiophagus.hannah     <- reportBestMatches(Ophiophagus.hannah.50hits)
#' best.hits.Pantherophis.guttatus  <- reportBestMatches(Pantherophis.guttatus.50hits)
#' 
#' ### Extract the subject sequences for the best matches
#' Crotalus.horridus.best.hits.seqs   <- get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.horridus,input.seqs=Crotalus.horridus.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Crotalus.horridus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' ### Same as next line:
#' Crotalus.horridus.best.hits.seqs   <- get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.horridus,input.seqs="/Users/alyssaleinweber/Documents/genomes/genomes_seqs/GCA_001625485.1_ASM162548v1_genomic.fna.gz",output.path="/Users/alyssaleinweber/Documents/REES_test_output/Crotalus.horridus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Crotalus.mitchellii.best.hits.seqs <- get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.mitchellii,input.seqs="/Users/alyssaleinweber/Documents/genomes/genomes_seqs/GCA_000737285.1_CrotMitch1.0_genomic.fna.gz",output.path="/Users/alyssaleinweber/Documents/REES_test_output/Crotalus.mitchellii_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Ophiophagus.hannah.best.hits.seqs  <- get.seqs.from.blastTable(input.blastTable=best.hits.Ophiophagus.hannah,input.seqs="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/516/915/GCA_000516915.1_OphHan1.0/GCA_000516915.1_OphHan1.0_genomic.fna.gz",output.path="/Users/alyssaleinweber/Documents/REES_test_output/Ophiophagus.hannah_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Pantherophis.best.hits.seqs        <- get.seqs.from.blastTable(input.blastTable=best.hits.Pantherophis.guttatus,input.seqs=Pantherophis.guttatus.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Pantherophis.guttatus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#'
#' param table.headers Logical indicating if the first line of input.blastTable is the column names. Argument ignored if input.blastTable is a data table, data frame, or character matrix.

#' @title Get UCEs From BLAST Table
#' 
#' Extracts the sequences defined in an input table (a blast table with NCBI format 6) from an input genome. The input genome should be the same genome used as the subject in the blast table.
#' If the query sequences in the blast table correspond to Ultraconserved Elements (UCEs), then the output is the UCEs of the species with input genome.
#' The quality of the output UCE dataset depends on the quality/completeness of the blast search.
#' 
#' @param genome.filepath Path to the input genome.
#' @param input.blastTable Path to the blast table containing best matches to query UCE sequences (from another species or individual) found in the input genome. This input table is the output of the reportBestMatches function. 
#' @param output.path Where to save the UCEs.
#' @return UCE sequences (or any extracted sequences of interest) written to a file in fasta format.
#' @export get.UCEs.from.blastTable
get.UCEs.from.blastTable <- function(genome.filepath,input.blastTable,output.path) {
	data.best           <- data.table::fread(file=input.blastTable,sep=",",header=T)
	query.subject.id    <- paste(as.character(data.best$qseqid),"_Subject=",as.character(data.best$sseqid),"_",as.character(data.best$sstart),"_",as.character(data.best$send),sep="")
	subject.id          <- gsub(".+_","",as.character(data.best$sseqid))
	subject.start       <- data.best$sstart
	subject.end         <- data.best$send
	### create an index of file 'foo.fasta'; this avoids having to actually copy or move the file to a new directory
	Rsamtools::indexFa(genome.filepath)
	fa <- Rsamtools::FaFile(genome.filepath)
	gr <- as(GenomeInfoDb::seqinfo(fa), "GRanges")
	
	scaff.matches.all   <- match(subject.id, names(gr))
	start.all           <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=min)
	end.all             <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=max)
	subranges           <- IRanges::IRanges(start=start.all ,end=end.all,names=subject.id)
	gsubranges          <- GenomicRanges::GRanges(seqnames=subject.id,ranges=subranges)
	UCE.scaff           <- Biostrings::getSeq(fa, gsubranges)
	names(UCE.scaff)    <- query.subject.id
	Biostrings::writeXStringSet(x = UCE.scaff, filepath=output.path, append=F, format="fasta")
}





