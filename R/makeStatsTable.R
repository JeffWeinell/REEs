#' @title Make Stats Table
#' 
#' This function aligns and then creates a table of stats calculated for each exon alignment. The output stats table is used to identify/select the REEs.
#' This function should be used instead of running the functions align.and.concatenate.best.exons and then alignment.stats
#' 
#' @param input.seqs Either a list of DNAStringSet objects with each species sequences, or a vector of character strings with paths to input sequences in fasta format.
#' @param species.names Vector of character strings containing names of species associated with each input.seqs.
#' @param reference.species Number indicating which species is the reference species.
#' @param input.gff One of the following: An object of class data table, data.frame, or character matrix, with columns identical to those used in Generic Feature Format (GFF) tables; or, a character string with local file path to either a GFF file or table with columns identical to those used in GFF files; or, a character string with URL path to a GFF file; or, NULL, in which case the gene.names column of the output table will be filled with NAs.
#' This table can either be the original GFF table or the one generated by the function filter.gff. Exon and protein isoform information may included in future versions of the output table, and this would require the original GFF table.
#' If not NULL, the GFF must correspond to the reference species.
#' @param subgroup Optional vector of character strings containing a subset of the "species" argument. Default is NULL. This parameter determines which species are used for calculating stats. If NULL, all species in the species argument are used.
#' @param output.path Where to save the stats table. If NULL, the stats table is saved to a temporary file.
#' @param alignments.out Directory where to save output alignment. Default is NULL (alignments not saved). If saved, alignments are in fasta format.
#' @param i.start First locus to start at. Default is 1.
#' @param i.stop Last locus to include. Default is NA, meaning to include all loci from i.start.
#' @param mafft.params Character string of mafft parameter settings. Default is "--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6". I don't suggest changing these.
#' @return Stats table (class data.table object), which is also saved to the value of output.path.
#' @export makeStatsTable
makeStatsTable <- function(input.seqs,species.names,reference.species,input.gff,output.path = NULL,alignments.out=NULL,subgroup=NULL,i.start=1,i.stop=NA,mafft.params="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6"){
	### Reorders the list of species such that the primary.species is first in the list
	species  <- c(species.names[reference.species],species.names[-reference.species])
	species  <- mgsub(c(" ","\\.","-"),c("_","_","_"),species)
	if(is.null(subgroup)){
		subgroup <- species
	}
	if(is.null(output.path)){
		output.path <- tempfile()
		delete <- T
	} else {
		delete <- F
	}
	### If alignments.out is not NULL, check if it exists and if not, create it.
	if(!is.null(alignments.out)){
		# dir.check.create(alignments.out)
		dir.create(alignments.out,recursive=T)
	}
	### A list of numbers indicating which of the species are also in the subgroup.
	is.subgroup     <- mgrep(query=subgroup,subject=species)
	if(is(input.seqs[[1]],"character")){
		### Puts the input.seqs in the same order as species.
		input.seqs   <- c(input.seqs[reference.species],input.seqs[-reference.species])
		### Read in the sequences for each species as a DNAStringSet, and hold the set of DNAStringSets in a list
		# Important not to name each DNAStringSet in the list
		species.seqs <- lapply(X=input.seqs,FUN=function(x){Biostrings::readDNAStringSet(x)})
	} else {
		if(is(input.seqs[[1]],"DNAStringSet")){
			##### Important not to name each DNAStringSet in the list
			### Puts the input.seqs in the same order as species.
			input.seqs        <- input.seqs[c(reference.species,c(1:length(input.seqs))[-reference.species])]
			### Set species.seqs equal to input.seqs in this case.
			species.seqs <- input.seqs
		}
		else {
			stop("input.seqs must be a list of DNAStringSet objects or a vector of character strings with filepaths to sequences in fasta format")
		}
	}
	### Extracts the QUERY contig accession name and range from each sequence name of each species.
	matches.list    <- lapply(X=species.seqs,FUN=function(x){ gsub("_Subject=.*","",names(x))})
	### This is the set of all loci found in at least one individual.
	all.loci        <- unique(unlist(matches.list))
	### List of numeric vectors, each vector containing numbers indicating which sequences are those in all.sequences
	matches.indices.all         <- lapply(X=matches.list,FUN=function(y){match(x=all.loci,table=y)})
	### Puts loci in the same order for each species. Doesnt filter based upon whether loci are shared. It's possible that species.seqs are already in the same order.
	species.seqs.all            <- mapply(FUN=function(X,Y){Z=Y[which(!is.na(Y))];X[Z]},X=species.seqs,Y=matches.indices.all,SIMPLIFY=F)
	### Create a vector holding contig names of all loci in input.seqs.
	CDS.locus.identifier.all    <- mgsub(c(":","-","CDS_","_Gene=.+"),c("_","_","",""),all.loci)
	### character string location of CDS within contig, with format "start_end"
	loci.ranges.all  <- mgsub(c(".*\\.1.",".*\\.1:","-","_Gene=.+"),c("","","_",""),all.loci)
	### length of each sequence (in species associated with the gff table) for the shared sequences
	loci.lengths <- (abs(as.numeric(gsub(".*_","",loci.ranges.all))-as.numeric(gsub("_.*","",loci.ranges.all)))+1)

	### Check if input.gff is NULL and if so set gene names to NA. Otherwise, read the GFF and extract gene names.
	if(is.null(input.gff)){
		gene.names.all <- rep(NA,length(all.loci))
	} else {
		if(is(input.gff,"data.table") | is(input.gff,"matrix")){
			annotationTable  <- as.data.frame(input.gff)
		}
		if(is(input.gff,"character")){
			### Reads in filtered GFF table, which is associated with the primary exome
#			annotationTable  <- data.table::fread(input=input.gff)
			annotationTable  <- as.data.frame(load.gff(input.gff))
		}
		##### Set column modes
		# Set which columns should be mode numeric
		numeric.columns <- which(colnames(annotationTable) %in% c("start","end"))
		# Set mode to numeric for those columns that should be numeric
		annotationTable[, numeric.columns] <- sapply(annotationTable[, numeric.columns], as.numeric)
		# Set which columns should be mode character
		character.columns <- which(!(colnames(annotationTable) %in% c("start","end")))
		# Set mode to "character" for the columns indexed in the character.columns vector
		annotationTable[, character.columns] <- sapply(annotationTable[, character.columns], as.character)
		## Vector of identifiers that should match those in CDS.locus.identifier.all
		annotationTable.identifier  <- paste0(unlist(annotationTable[,1]),"_",unlist(annotationTable[,"start"]),"_", unlist(annotationTable[,"end"]))
		## Extract gene name from the last column of the gff table
#		gene.names.temp             <- mgsub(c(".*;gene=",";.*"),c("",""),unlist(annotationTable[,9]))
#		gene.names.temp2            <- gsub(";.*","",gene.names.temp)
#		match.identifier.all        <- match(CDS.locus.identifier.all,annotationTable.identifier)
		### Potential alternative to the previous line, because the match function only uses first feature annotated for each CDS.
		match.identifier.all        <- lapply(X=paste0("^",CDS.locus.identifier.all,"$"),FUN=grep, annotationTable.identifier)
#		gene.names.all              <- unlist(lapply(X=match.identifier.all,FUN=function(x){paste0("gene=",paste0(unique(gene.names.temp[x]),collapse=","))}))
		# range.annotations.all is a list of lists of character vectors, each of which contains the annotations (column 9) for one feature of a genomic region in CDS.locus.identifier.all
		range.annotations.all       <- lapply(X=match.identifier.all,FUN=function(x){strsplit(annotationTable[x,"attributes"],split=";")})
		# range.annotations.all.compressed is a list of character vectors, each of which contains all annotations (column 9) for all features with of a genomic region in CDS.locus.identifier.all
		range.annotations.all.compressed <- lapply(range.annotations.all,FUN=unlist)
		gene.names.all <- unlist(lapply(X=range.annotations.all.compressed,FUN=function(x){paste0("gene=",paste(unique(gsub("^gene=","",x[grep("^gene=",x)])),collapse=","))}))
	}

	if(is.na(i.start) | i.start > length(all.loci)){
		### Default i.start, ie which exons to start the loop at
		i.start <- 1
	}
	if(is.na(i.stop) | i.stop > length(all.loci)){
		### Default i.stop, ie which exons to stop the loop at
		i.stop  <- length(all.loci)
	}
	# List of unaligned DNAStringSet objects. Loci do not need to be shared.
	index.matrix.all            <- do.call(rbind, matches.indices.all)
	# The function defining 'temp.seqs' object will only work if names(species.seqs.all) is NULL, so the next line sets names to NULL if they arent already NULL
	if(!is.null(names(species.seqs.all))){
		names(species.seqs.all) <- NULL
	}
	# This will take a few minutes to run
	temp.seqs               <- lapply(c(1:ncol(index.matrix.all)),function(input){do.call(c,mapply(FUN=function(A,B){C=B[which(!is.na(B))];A[C]},A=species.seqs.all,B=lapply(X=index.matrix.all[,input],FUN=function(x){x})))})
	# Matrix or list of species names to use for each sequence in each alignment.
	names.temp.seqs         <- apply(X=index.matrix.all,MARGIN=2,FUN=function(input){species[!is.na(input)]})
	# Now setting the species name to each sequence in each DNAStringSet
	if(is(names.temp.seqs,"list")){
		for(i in 1:length(temp.seqs)){
			names(temp.seqs[[i]]) <- names.temp.seqs[[i]]
		}
	}
	if(is(names.temp.seqs,"matrix")){
		for(i in 1:length(temp.seqs)){
			names(temp.seqs[[i]]) <- names.temp.seqs[,i]
		}
	}
	### Create a matrix to hold alignment stats for all loci
	tempMatrix           <- matrix(data=0, nrow=length(1:length(all.loci)), ncol=12+length(species))
	colnames(tempMatrix) <- c(paste0(species[1],".locus"), "num.Species","CountCover","absolutePIS","percentPIS","mean.pident", paste0("pident.",species), "gene.name",paste0("locus.length.",species[1]),"mean.variable.sites","min.pident.all","min.pident.subgroup","alignment.width")
	
	# Coerce tempMatrix to a data frame?
	# tempMatrix <- as.data.frame(tempMatrix)
	### Set column modes. All columns should be numeric except the first column and the "gene.name", which should be "character" mode.
	# Find which column is the gene.name column
	which.is.gene.name.column <- which(colnames(tempMatrix)=="gene.name")
	# Columns that should be character mode
	character.columns <- c(1,which.is.gene.name.column)
	# Set the mode to numeric for those columns that should be character
	tempMatrix[, character.columns] <- sapply(tempMatrix[, character.columns], as.character)
	# Columns that should be numeric mode
	numeric.columns   <- setdiff(1:ncol(tempMatrix),character.columns)
	# Set the mode to numeric for those columns that should be numeric
	tempMatrix[, numeric.columns] <- sapply(tempMatrix[, numeric.columns], as.numeric)
	### Set column 1 to 
	tempMatrix[,1]             <- all.loci
	tempMatrix[,"gene.name"]   <- gsub("gene=","",gene.names.all)
	tempMatrix[,"num.Species"] <- sapply(temp.seqs,length)
	### Will be the same for all species because intersect.all function was used.
	### This currently only considers tempMatrix and not tempMatrix
	for(i in c(i.start:i.stop)) {
		### Make an alignment for ith locus if more than one sequence in ith DNAStringSet
		if (tempMatrix[i,"num.Species"]>1){
			### aligns the ith locus of each species
#			alignment.temp              <- REEs::mafft(temp.seqs[[i]],param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
			alignment.temp              <- REEs::mafft(temp.seqs[[i]],param=mafft.params)
			### Number of sites in alignment
			alignment.width             <- width(alignment.temp)[1]
			### pairwise matrix of absolute genetic distance
			distances                   <- 100*as.matrix(ape::dist.dna(ape::as.DNAbin(alignment.temp),model="raw",pairwise.deletion=T))
			### Dont use Biostrings::stringDist method for calculating genetic distances because it doesnt allow for pairwise deletions. If it did then the "hamming" method would be the same as ape's "raw" method.
			# Biostrings::stringDist pairwise applies Biostrings::neditStartingAt.
			# distances <- (100*(Biostrings::stringDist(alignment.temp)/alignment.width))
			### pairwise matrix of percent genetic identity
			pident                      <- round(100-distances,digits=2)
			mean.pident                 <- round(mean(pident[1,-1]),digits=2)
			### mean percent identity of the primary species to each of the other species.
			tempMatrix[i,"mean.pident"] <- mean.pident
			### percent identity of the primary species to each species in the alignment (including itself, with should be 100% idendical).
			tempMatrix[i,(7:(length(species)+6))[species %in% names(alignment.temp)]] <- pident[1,]
			if(any(!species %in% names(alignment.temp))){
				tempMatrix[i,(7:(length(species)+6))[!species %in% names(alignment.temp)]] <- NA
			}
		}
		if (tempMatrix[i,"num.Species"] > 3){
			# Counts the number of sites of ith locus with at least 4 individuals with non-missing data
			count.cover <- width(filter.alignment(alignment=alignment.temp,treat.ambiguous.as.missing=T,min.allele.freqs.dna=c(4,0,0,0)))[1]
		} else {
			count.cover <- 0
		}
		# If no sites have ≥4 individuals with non-missing data,then the number and percent of parsimony informative sites is zero for ith locus.
		# Otherwise, calculates the number of parsimony informative sites (using pis function) for ith locus
		if(is.na(count.cover) | is.null(count.cover) | count.cover == 0){
			count.cover <- 0
			absolutePIS <- 0
			percentPIS  <- 0
		} else {
			absolutePIS <- REEs::pis(ape::as.DNAbin(Biostrings::DNAMultipleAlignment(alignment.temp)),what="absolute")
			percentPIS  <- round(absolutePIS/count.cover,digits=4)*100
		}
		tempMatrix[i,"CountCover"]      <- count.cover  ### number of characters at locus with at least 4 individuals represented
		tempMatrix[i,"absolutePIS"]     <- absolutePIS  ### number of parsimony informative sites at locus
		tempMatrix[i,"percentPIS"]      <- percentPIS   ### percent of characters with at least 4 individuals represented that are parsimony informative
		tempMatrix[i,"alignment.width"] <- alignment.width
		### Calculates the mean number of variable sites; this value isnt used later
		mean.var.sites      <- round(((100-as.numeric(mean.pident))/100)*as.numeric(loci.lengths[i]),digits=2)
		### Minimum pident of any species to the primary species for ith locus (excluding species without data)
		min.pident.all      <- round(min(as.numeric(tempMatrix[i,(7:(length(species)+6))[species %in% names(alignment.temp)]])),digits=2)
		### Minimum pident of any species in the subgroup to the primary species for ith locus.
		if(all(subgroup %in% species)){
			min.pident.subgroup <- min.pident.all
		} else {
			min.pident.subgroup <- min(as.numeric(tempMatrix[(is.subgroup+6)]))
		}
		tempMatrix[i,paste0("locus.length.",species[1])] <- loci.lengths[i]
		tempMatrix[i,"mean.variable.sites"] <- mean.var.sites
		tempMatrix[i,"min.pident.all"]      <- min.pident.all
		tempMatrix[i,"min.pident.subgroup"] <- min.pident.subgroup
		alignment.out <- Biostrings::DNAStringSet(alignment.temp)
		if(!is.null(alignments.out)){
			Biostrings::writeXStringSet(x=alignment.out,filepath=paste0(alignments.out,"/locus",i,"_aligned.fas"))
		}
	}
	### writes the column names and tempMatrix for exon i = 1
	write.table(x=tempMatrix,file=output.path,sep="\t",col.names=T,row.names=F,append=F)
	if(delete){
		file.remove(output.path)
	}
	result <- data.table::as.data.table(tempMatrix)
	result
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
#' ### Use get.seqs.from.gff to extract the sequences for the loci in the filtered GFF
#' Thamnophis.sirtalis_genome.path <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz"
#' Thamnophis.sirtalis_exome       <- get.seqs.from.gff(input.seqs=Thamnophis.sirtalis_genome.path,input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp)
#' 
#' ### Use tblastx to return up to 50 matches for each of query sequence in a subject sequence database
#' # Define the query sequences. In this example, query sequences are the first sequences in the Thamnophis sirtalis exome DNAStringSet object.
#' test.query    <- Thamnophis.sirtalis_exome[1:2]
##'### Or read from fasta file
#' test.query    <- readDNAStringSet("/Users/alyssaleinweber/Documents/REES_test_output/Thamnophis.sirtalis_twoExons_testQuery.fas")
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
#' Gekko.japonicus.50hits              <- blast(method="tblastx",subject=Gekko.japonicus.genome_url,query=test.query,table.out="/Users/alyssaleinweber/Documents/REES_test_output/Gekko.japonicus_TwoExons.testQuery.tblastx.50hits.txt")
#' Pagona.vitticeps.50hits             <- blast(method="tblastx",subject=Pagona.vitticeps.genome_url,query=test.query,table.out="/Users/alyssaleinweber/Documents/REES_test_output/Pagona.vitticeps_TwoExons.testQuery.tblastx.50hits.txt")
#' Crotalus.horridus.50hits            <- blast(method="tblastx",subject=Crotalus.horridus.genome_url,query=test.query)
#' Crotalus.mitchellii.50hits          <- blast(method="tblastx",subject=Crotalus.mitchellii.genome_url,query=test.query)
#' Ophiophagus.hannah.50hits           <- blast(method="tblastx",subject=Ophiophagus.hannah.genome_url,query=test.query)
#' Pantherophis.guttatus.50hits        <- blast(method="tblastx",subject=Pantherophis.guttatus.genome_url,query=test.query)
#' Protobothrops.mucrosquamatus.50hits <- blast(method="tblastx",subject=Protobothrops.mucrosquamatus.genome_url,query=test.query)
#' Python.bivittatus.50hits            <- blast(method="tblastx",subject=Python.bivittatus.genome_url,query=test.query)
#' Viperus.berus.50hits                <- blast(method="tblastx",subject=Viperus.berus.genome_url,query=test.query)
#' 
#' ### Filter results to include only the best match for each query sequence
#' best.hits.Anolis.carolinensis           <- reportBestMatches(Anolis.carolinensis.50hits)
#' best.hits.Gekko.japonicus               <- reportBestMatches(Gekko.japonicus.50hits,output.table.path="/Users/alyssaleinweber/Documents/REES_test_output/Gekko.japonicus_TwoExons.testQuery.tblastx.best.hits.txt")
#' best.hits.Pagona.vitticeps              <- reportBestMatches(Pagona.vitticeps.50hits,output.table.path="/Users/alyssaleinweber/Documents/REES_test_output/Pagona.vitticeps_TwoExons.testQuery.tblastx.best.hits.txt")
#' best.hits.Crotalus.horridus             <- reportBestMatches(Crotalus.horridus.50hits)
#' best.hits.Crotalus.mitchellii           <- reportBestMatches(Crotalus.mitchellii.50hits)
#' best.hits.Ophiophagus.hannah            <- reportBestMatches(Ophiophagus.hannah.50hits)
#' best.hits.Pantherophis.guttatus         <- reportBestMatches(Pantherophis.guttatus.50hits)
#' best.hits.Protobothrops.mucrosquamatus  <- reportBestMatches(Protobothrops.mucrosquamatus.50hits)
#' best.hits.Python.bivittatus             <- reportBestMatches(Python.bivittatus.50hits)
#' best.hits.Viperus.berus                 <- reportBestMatches(Viperus.berus.50hits)
#' 
#' ### Extract the subject sequences for the best matches
#' Crotalus.horridus.best.hits.seqs <- get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.horridus,input.seqs=Crotalus.horridus.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Crotalus.horridus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' ### Same as next line:
#' Anolis.carolinensis.best.hits.seqs          <- get.seqs.from.blastTable(input.blastTable=best.hits.Anolis.carolinensis,input.seqs=Anolis.carolinensis.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Anolis.carolinensis_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Gekko.japonicus.best.hits.seqs              <- get.seqs.from.blastTable(input.blastTable=best.hits.Gekko.japonicus,input.seqs=Gekko.japonicus.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Gekko.japonicus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Pagona.vitticeps.best.hits.seqs             <- get.seqs.from.blastTable(input.blastTable=best.hits.Pagona.vitticeps,input.seqs=Pagona.vitticeps.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Pagona.vitticeps_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Crotalus.horridus.best.hits.seqs            <- get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.horridus,input.seqs="/Users/alyssaleinweber/Documents/genomes/genomes_seqs/GCA_001625485.1_ASM162548v1_genomic.fna.gz",output.path="./Crotalus.horridus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Crotalus.mitchellii.best.hits.seqs          <- get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.mitchellii,input.seqs="/Users/alyssaleinweber/Documents/genomes/genomes_seqs/GCA_000737285.1_CrotMitch1.0_genomic.fna.gz",output.path="./Crotalus.mitchellii_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Ophiophagus.hannah.best.hits.seqs           <- get.seqs.from.blastTable(input.blastTable=best.hits.Ophiophagus.hannah,input.seqs="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/516/915/GCA_000516915.1_OphHan1.0/GCA_000516915.1_OphHan1.0_genomic.fna.gz",output.path="./Ophiophagus.hannah_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Pantherophis.guttatus.best.hits.seqs        <- get.seqs.from.blastTable(input.blastTable=best.hits.Pantherophis.guttatus,input.seqs=Pantherophis.guttatus.genome_url,output.path="./Pantherophis.guttatus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Protobothrops.mucrosquamatus.best.hits.seqs <- get.seqs.from.blastTable(input.blastTable=best.hits.Protobothrops.mucrosquamatus,input.seqs=Protobothrops.mucrosquamatus_genome.url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Protobothrops.mucrosquamatus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Python.bivittatus.best.hits.seqs            <- get.seqs.from.blastTable(input.blastTable=best.hits.Python.bivittatus,input.seqs=Python.bivittatus.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Python.bivittatus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' Viperus.berus.best.hits.seqs                <- get.seqs.from.blastTable(input.blastTable=best.hits.Viperus.berus,input.seqs=Viperus.berus.genome_url,output.path="/Users/alyssaleinweber/Documents/REES_test_output/Viperus.berus_TwoExons.testQuery.tblastx.best.hits_seqs.fas")
#' 
#' ### Align homologous sequences and make a stats table to summarizing variation in each alignment (one row per aligned locus).
#' input.seqs.paths <- .... ### paths to the "best.hits.seqs" files
#' stats.table      <- makeStatsTable(input.seqs=input.seqs.paths,species=species.temp,input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp,output.path=table.out,alignments.out=alignments.dir,reference.species=10)

#' @title Align and Concatenate REEs
#' 
#' Takes as input a set of exomes and an input table indicating which loci are the REEs, performs multiple sequence alignment (using mafft) for each REE, and outputs:
#' (1) An alignments for each REE locus,
#' (2) The concatenated alignment of REE loci,
#' (3) The partition file associated with the concatenated locus alignment.
#'
#' This function is superceded by makeStatsTable.
#'
#' @param input.seqs Paths to unaligned input sequences, one file per species. These files should be the output files generated by the function get.seqs.from.BlastTable()
#' @param statsTable.path Path to the stats table generated by the function makeExomeStatsTable; Example: "stats_data_FastestExonPerGene_best.txt"
#' @param output.dir Where to save the (1) the alignments for each locus, (2) the concatenated locus alignment, and (3) the partition file associated with the concatenated locus alignment. Example: "Exomes_TempFolder_3Nov2019/"
#' @param species Species names associated with individuals in the alignments.
#' @param is.primary.exome Number indicating which exome is the primary exome, i.e., the one that other exomes are aligned to and the one that divergence statistics are calculated against.
#' @param i.start First locus to start at. Default is 1.
#' @param i.stop Last locus to include. Default is NA, meaning to include all loci from i.start.
#' @return Separate DNA alignments for each locus in the input stats table, a concatenated-locus alignment, and a partition file for the concatenated locus alignment.
#' @export align.and.concatenate.best.exons
align.and.concatenate.best.exons <- function(input.seqs,statsTable.path,output.dir,species,is.primary.exome,i.start = 1,i.stop = NA){
	statsTable           <- data.table::fread(statsTable.path,sep=",",header=T)             ### reads in the info on fastest exon per gene dataset with total length <1.2Mbp
	species              <- c(species[is.primary.exome],species[-is.primary.exome])         ### reorders the list of species such that the primary.species is first in the list
	input.seqs     <- input.seqs[mgrep(species,input.seqs)]               ### puts the input.seqs in the same order as species
	exome.names          <- paste("exome",c(1:length(input.seqs)),sep="")
	
	for(i in 1:length(exome.names)){                                                       #|reads in exomes
		assign(x=exome.names[i],value=Biostrings::readDNAStringSet(filepath=input.seqs[i]))      #|and assigns them
	}                                                                                      #|to object names

	CDS.names        <- names(exome1)
	unique.CDS.names <- unique(CDS.names)
	to.keep          <- match(unique.CDS.names, CDS.names)
	CDS.names        <- unique.CDS.names
	exome1           <- exome1[to.keep]

	locus.names <- paste(unlist(statsTable[,1]),unlist(statsTable$gene.name),sep="_")
	CDS.names   <- unlist(statsTable[,1])

	if(is.na(i.start) | i.start > length(CDS.names)){
		i.start <- 1                  ### default i.start, ie which exons to start the loop at
	}
	if(is.na(i.stop) | i.stop > length(CDS.names)){
		i.stop  <- length(CDS.names)   ### default i.stop, ie which exons to stop the loop at
	}
	
	# dir.check.create(paste(output.dir,"exonAlignments",sep=""))
	dir.create(paste(output.dir,"exonAlignments",sep=""),recursive=T)
	
	for(i in i.start:i.stop){
	
		temp.exon.names <- paste("temp.exon.exome",c(1:length(species)),sep="")
		for(k in 1:length(temp.exon.names)){
			assign(x=temp.exon.names[k],value = get(exome.names[k])[grep(CDS.names[i],names(get(exome.names[k])))])
		}
		temp.seqs <- Biostrings::DNAStringSet(get(temp.exon.names[1]))
		for(z in 2:length(temp.exon.names)){
			temp.seqs <- c(temp.seqs,Biostrings::DNAStringSet(get(temp.exon.names[z])))
		}
		### Aligns the ith exon of each species
		alignment           <- REEs::mafft(temp.seqs,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
		Biostrings::writeXStringSet(x=alignment, filepath=paste(output.dir,"exonAlignments/",locus.names[i],".fas",sep=""), append = F)
		
		if(i == 1){
			alignment2 <- alignment            ### concatenated alignment for i = 1
			start.temp <- 1                    ### site where ith partition starts for i = 1
			end.temp   <- width(alignment)[1]  ### site where ith partition ends for i = 1
			app=FALSE                          ### states that we do not want to the append partition info to next line of partition file for the first loop
		}
		if(i != 1){
			start.temp <- width(alignment2)[1]+1                  ### width of the (i-1)th alignment2 + 1
			alignment2 <- Biostrings::xscat(alignment,alignment2) ### concatenates the ith alignment with the (i-1)th alignment2
			end.temp   <- width(alignment2)[1]                    ### width of the ith alignment2
			app=TRUE                                              ### states that we will want to the append partition info to next line of partition file
		}
		temp.write <- paste("DNA, ",locus.names[i]," = ",start.temp,"-",end.temp,sep="")
		temp.write <- gsub("Subject=","Subject_",temp.write)
		temp.write <- gsub("gene=","gene_",temp.write)
		write(temp.write, paste(output.dir,"fastestExonPerGene_BestExons_partitionFile.txt",sep="") ,append=app)
		if(i%%10==0){print(paste(i,"of",length(CDS.names),"complete",Sys.time(),sep=" "))}
	}
	names(alignment2) <- species
	Biostrings::writeXStringSet(alignment2, paste(output.dir,"FastestExonPerGene_BestExons_aligned_concatenated.txt",sep=""), append = F) #writes the aligned and concatenated target loci to file
}

#' @title Align Shared Loci
#' 
#' Takes as input a set of fasta sequence files, one per species (or individual). Perfoms multiple sequence alignment (using MAFFT) for loci that are shared among all species (determined by sequence names).
#' Returns a list of DNAStringSet objects, each of which holds aligned sequences for a locus. Optionally writes alignments in fasta format (one alignment per file).
#'
#' @param input.seqs Character vector of filepaths to the fasta sequence files. Each input file holds sequence data for one species (or for one individual if multiple individuals per species). Fasta headers (sequence names) must have a consistent, underscore-delimited format, with locus name being one of the entries (Default is for locus name to be held between the second underscore and either the end of the header or another underscore). Example header that I used for UCEs: "Genus_species_uce-XXXX", where XXXX is a number unique to a particular UCE homolog.
#' @param indv Character vector of names of individuals (or species if one individual per species) associated with input sequences. Order of names in indv must correspond to order of files in input.seqs, and each name in indv needs to occur in its corresponding filename in input.seqs.
#' @param reference.indv Number indicating which individual in indv is the reference (aka primary) individual for sequence alignment. Default is 1.
#' @param seqname.str.delim Character delimiter to use for parsing fasta headers of input.seqs. Default is "_". After parsing headers by this delimiter, one of the entries (set by the seqname.str.loc parameter) holds the locus name of the sequence.
#' @param seqname.str.loc A number specifying which part of fasta sequence header lines (after parsing by seqname.str.delim) contains the locus name. Default is 3. For example: ">Something_Anything_Locus1" or "Text_MoreText_Locus1_MoreInfo" would both be appropriate with default seqname.str.loc value of 3.
#' @param output.dir Where to save the output fasta alignment files. Default is NULL.
#' @param i.start First locus to start at. Default is 1. For now only works if set to 1.
#' @param i.stop Last locus to include. Default is NA, meaning to include all loci from i.start.
#' @return Separate DNA alignments for each locus in the input stats table, a concatenated-locus alignment, and a partition file for the concatenated locus alignment.
#' @export align.shared.loci
align.shared.loci <- function(input.seqs,indv,reference.indv=1,seqname.str.delim="_",seqname.str.loc=3,output.dir=NULL,i.start = 1,i.stop = NA){
	if(i.start!=1){
		stop("in the current version of this function, i.start must be set to 1")
	}
	### Reorders the list of species such that the reference individual is first in the list
	indv         <- c(indv[reference.indv],indv[-reference.indv])
	### Puts the input.seqs in the same order as species.
	input.seqs   <- input.seqs[mgrep(indv,input.seqs)]
	### Reads in the sequences for each species as a DNAStringSet, and holds the set of DNAStringSets in a list
	species.seqs <- lapply(X=input.seqs,FUN=function(x){Biostrings::readDNAStringSet(x)})
	### List of vectors, with the ith vector holding the names of loci (extracted from sequence names) of the ith individual in input.seqs
	loci.list  <- lapply(X=species.seqs,FUN=function(x){vapply(strsplit(names(x),seqname.str.delim), `[`, seqname.str.loc, FUN.VALUE=character(1))})
	### Table indicating number of species that have each locus
	loci.table    <- table(unlist(loci.list))
	### Vector of UCE names for UCEs present in all species in input.seqs
	shared.loci   <- names(loci.table)[which(loci.table == length(indv))]
	if(length(shared.loci)==0){
		stop("Zero loci shared by all individuals")
	}
	### A list of vectors. The jth element of the ith vector indicates which sequence of the ith species is the jth locus in shared.loci
	matches.indices.shared <- lapply(X=loci.list,FUN=function(y){match(x=shared.loci,table=y)})
	### species.seqs filtered to only include shared loci, and sorted such that sequences are in the same order for each species.
	species.seqs.shared    <- mapply(FUN=function(X,Y){Z=Y[which(!is.na(Y))];X[Z]},X=species.seqs,Y=matches.indices.shared,SIMPLIFY=F)
	if(is.na(i.stop)){
		# i.stop  <- length(UCEs.shared)
		i.stop  <- length(shared.loci)
	}
	# List of unaligned DNAStringSet objects. Loci do not need to be shared.
	temp.seqs <- sapply(X=c(i.start:i.stop),FUN=function(z){seqset=do.call(c,mapply(FUN=function(X,Y){X[Y]},X=species.seqs.shared,Y=rep(list(z),length(indv)),SIMPLIFY=F));names(seqset)=species;seqset})
	### Create an empty vector that will be filled with the alignments
	result.list <- list(); length(result.list) <- length(i.start:i.stop)
	### Align each locus!
	for(i in i.start:i.stop){
		### Aligns the ith locus of each species
		alignment.temp           <- REEs::mafft(temp.seqs[[i]],param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
		if(!is.null(output.dir)){
			### Define output path for the alignment of the ith locus.
			output.path.temp <- gsub("//","/",paste0(output.dir,"/",shared.loci[i],".fas"))
			### Save alignment of ith locus
			Biostrings::writeXStringSet(alignment.temp, output.path.temp, append = F)
		}
		result.list[[i]]        <- alignment.temp
		names(result.list)[[i]] <- shared.loci[i]
	}
	#result.list      <- result.list
#	alignment.widths <- unlist(lapply(X=result.list,FUN=function(x){width(x)[1]}))
#	loci.end         <- rollSum(alignment.widths)
#	loci.start       <- (c(0,loci.end)+1)[1:length(alignment.widths)]
#	ranges.mat       <- cbind(names(result.list),loci.start,loci.end)
	### Concatenate loci. Names will be lost temporarily.
#	result.concat <- do.call(xscat,result.list)
	### Reset sequence names for the concatenated loci alignment as the species names.
	#names(result.concat) <- names(result.list[[1]])
	result.list
}


#' @title Get Alignment Stats Table
#' 
#' Similar to the function makeExomeStatsTable, except this function takes as input the set of exon alignments rather than the unaligned exomes. The output is an abbreviated version of the stats table produce by the function makeExomeStatsTable.
#' This function is superceded by makeStatsTable have been abandoned. Note that pis.new is used here rather than the pis function that was used in makeExomeStatsTable function.
#' 
#' @param align.dir Directory where exon alignments are located.
#' @param outfile Path where the output stats table should be saved.
#' @param seqs.format File format of input exon alignments.
#' @return Abbreviated exome stats table writen to a file. The stats table is also returned as a character matrix.
#' @export alignment.stats
alignment.stats <- function(align.dir,outfile,seqs.format="phylip"){
	ext         <- c("phy","fa"); names(ext) <- c("phylip","fasta")
	ext.temp    <- as.character(ext[seqs.format])
	### List of the names of alignment files
	files       <- list.files(path=align.dir, pattern=paste("*.",ext.temp,sep=""),full.names=T)
	### List of the names of alignment files
	short.names <- list.files(path=align.dir, pattern=paste("*.",ext.temp,sep=""))
	stats.table <- matrix(data="0",nrow=length(files),ncol=10)
	colnames(stats.table) <- c("locus","n_individuals","n_characters","percent_gaps","fraction_sites_pars_inform","mean_pDistance","all_some_overlap","stat7","stat8","stat9","stat10")[1:ncol(stats.table)]
	for(i in 1:length(files)){
		
		multipleAlignment.temp <- readDNAMultipleAlignment(file = files[i], format = seqs.format)
		### DNAbin version of the alignment
		DNAbin.temp    <- Biostrings::as.DNAbin(multipleAlignment.temp)
		alignment.temp <- Biostrings::unmasked(multipleAlignment.temp)
		
		distances        <- ape::dist.dna(DNAbin.temp,model="raw",pairwise.deletion=T)
		mean.distance    <- round(mean(distances[which(distances!="NaN")]),digits=3)
		all.some.overlap <- all(distances!="NaN")
		name.temp        <- gsub(pattern=".phy","",x=short.names[i])
		stats.table[i,1] = name.temp
		stats.table[i,2] = as.numeric(length(alignment.temp))
		stats.table[i,3] = as.numeric(width(alignment.temp)[1])
		### Percent of sites that are informative
		stats.table[i,5] = pis.new(alignment.temp,abs=F)
		stats.table[i,6] = mean.distance
		stats.table[i,7] = all.some.overlap
		
		if(i==1){
			write(x=colnames(stats.table),ncolumns=ncol(stats.table),file=outfile,append=T)
		}
		write(x=stats.table[i,],ncolumns=ncol(stats.table),file=outfile,append=T)
	}
	stats.table
}
