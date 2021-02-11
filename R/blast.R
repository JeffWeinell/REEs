#' @title Run blast from R 
#' 
#' Wrapper for running BLAST from R 
#' Requires a local copy of BLAST, which can be dowloaded from here: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#' If blasting against genomes, this function is best implemented on a cluster/cloud environment, rather than locally, unless you have a super-computer...
#' 
#' @param blast.path Either "auto" (the default), or the full path to the directory containing blast executables (usually the bin subfolder of a folder named "ncbi-blast-#.##.#+", where # are numbers indicating the version).
#' "auto" indicates that blast was installed to the REEs package folder using the function blast.install(). Future versions of this function will seach for executables on your system path, but this is not currently implemented.
#' @param method One of the following: "blastn", "blastp", "blastx", "deltablast", "psiblast", "rpsblast", "rpstblastn", "tblastn", "tblastx"
#' @param subject One of the following: 
#' Character string indicating path to a local copy of a genome.
#' Character string indicating URL to fasta-formatted sequences.
#' DNAStringSet object. If a DNAStringSet object is used, then a temporary fasta file is created to hold the sequences to search within.
#' If DNAStringSet or URL character string, a local copy will be saved as a temporary file.
#' @param query Path to the fasta file containing the query sequences to find matches for.  If a DNAStringSet object is used, then a temporary fasta file is created to hold the query sequences.
#' @param table.out Where to save the output table. If NULL, the output is written to a temporary file that is read into R.
#' @param eval Number to use for the expect value (e-value); default is 1e-5.
#' @param output.format Integer indicating which format to use for the output table of matches (default 6).
#' @param max.targets.per.query Integer indicating the maximum number of targets (subject contigs) containing a match per query sequence. This sets the BLAST's max_target_seqs argument. Importantly, this argument does not control the total number of matches per query, because multiple matches may occer on the same subject contig.
#' @param max.matches.per.target Integer indicating the maximum number of matches per target (subject contig) per query. This sets BLAST's max_hsps argument. To determine the maximum total number of matches per query, multiply max.targets.per.query by max.matches.per.target.
#' @param parallel.groups  Number of groups to split the query sequences into for running in parallel. This is experimental.
#' @param num.threads Either an integer indicating how many threads to use, or "max" (default), in which case num.threads is set to the number of cores available.
#' @param other.args A character string of the form "-argument1 value1 -argument2 value2" with additional arguments and values to pass to BLAST. See BLAST manual for definitions of available arguments. Default is NULL.
#' @return Table of matches.
#' @export blast
blast <- function(blast.path="auto",method,subject,query,table.out=NULL,eval=1e-5,output.format=6,max.targets.per.query=10,max.matches.per.target=10,parallel.groups=NULL,num.threads="max",other.args=NULL){
	#### Prepare the path to the executables
	if(blast.path=="auto"){
		REEs.blast.dir   <- paste0(find.package("REEs"),"/blast-mafft/blast")
		blast.dir.path   <- list.dirs(REEs.blast.dir)[grep("bin$",list.dirs(REEs.blast.dir))]
	} else {
		### In the future, include another if statement to check if any directories on the path contain blast executables.
		blast.dir.path <- blast.path
	}
	### path to the executable of whichever blast method is to be implemented
	blast.exe.path       <- paste(blast.dir.path,method,sep="/")
	### path to the executable makeblastdb
	makeblastdb.exe.path <- paste(blast.dir.path,"makeblastdb",sep="/")
	#### Verify that blast.exe.path and makeblastdb.exe.path are executable. Stop and warn if not.
	test.blast.exe       <- check.if.executable(exe.path=blast.exe.path)
	test.makeblastdb.exe <- check.if.executable(exe.path=makeblastdb.exe.path)
	if(test.blast.exe!=0){
		stop(paste0("'",blast.exe.path," is not executable. Aborting. Run blast.install() with argument defaults to install BLAST to REEs package.'"))
	}
	if(test.makeblastdb.exe!=0){
		stop(paste0("'",makeblastdb.exe.path," is not executable. Aborting.'"))
	}
	#### Make a temporary fasta file holding the subject sequences if the value of the "subject" argument is a DNAStringSet or URL string.
	if("DNAStringSet" %in% class(subject)){
		subject.path   <- tempfile()
		names(subject) <- gsub(" .+","",names(subject))
		Biostrings::writeXStringSet(x = subject, filepath=subject.path, append=F, format="fasta")
		delete.subject <- T
	} else {
		if(file.exists(subject)){
			subject.path     <- subject
			subject.obj.temp <- Biostrings::readDNAStringSet(subject.path)
			names(subject.obj.temp) <- gsub(" .+","",names(subject.obj.temp))
			subject.path <- tempfile()
			Biostrings::writeXStringSet(x = subject.obj.temp, filepath=subject.path, append=F, format="fasta")
			delete.subject <- T
			rm(subject.obj.temp)
		} else {
			subject.path <- tempfile()
			# sets time limit for downloading files to 1000 seconds
			options(timeout=1000)
			conn <- utils::download.file(url=subject, destfile=subject.path)
			subject.obj.temp        <- Biostrings::readDNAStringSet(subject.path)
			names(subject.obj.temp) <- gsub(" .+","",names(subject.obj.temp))
			Biostrings::writeXStringSet(x = subject.obj.temp, filepath=subject.path, append=F, format="fasta")
			delete.subject <- T
			rm(subject.obj.temp)
		}
	}
	#### Make a temporary fasta file holding the query sequences if the value of the query argument is a DNAStringSet or URL string.
	if("DNAStringSet" %in% class(query)){
		query.path   <- tempfile()
		names(query) <- mgsub(c(" ",","),c("_","_"),names(query))
		names(query) <- substring(names(query),first=1,last=50)
		if(is.null(parallel.groups)){
			Biostrings::writeXStringSet(x = query, filepath=query.path, append=F, format="fasta")
			delete.query <- T
		} else {
			temp.file <- list()
			length(temp.file) <- parallel.groups
			out.files.temp    <- temp.file
			for(i in 1:parallel.groups){
				temp.file[[i]]      <- tempfile()
				out.files.temp[[i]] <- tempfile()
			}
			temp.file      <- unlist(temp.file)
			out.files.temp <- unlist(out.files.temp)
			groups         <- sort(sample(1:parallel.groups, size=length(query),replace=T))
			query.groups   <- lapply(X=c(1:parallel.groups),FUN=function(x){query[which(groups==x)]})
			for(i in 1:parallel.groups){
				Biostrings::writeXStringSet(x = query.groups[[i]], filepath=temp.file[[i]], append=F, format="fasta")
			}
		}
	} else {
		if(file.exists(query)){
			query.path     <- query
			query.obj.temp <- Biostrings::readDNAStringSet(query.path)
			if(any(nchar(unlist(names(query.obj.temp)))>50) | any(!is.na(stringr::str_locate(unlist(names(query.obj.temp))," ")))){
				query.path <- tempfile()
				names(query.obj.temp) <- mgsub(c(" ",","),c("_","_"),names(query.obj.temp))
				names(query.obj.temp) <- substring(names(query.obj.temp),first=1,last=50)
				if(is.null(parallel.groups)){
					Biostrings::writeXStringSet(x = query.obj.temp, filepath=query.path, append=F, format="fasta")
					delete.query <- T
				} else {
					temp.file <- list()
					length(temp.file) <- parallel.groups
					out.files.temp <- temp.file
					for(i in 1:parallel.groups){
						temp.file[[i]]      <- tempfile()
						out.files.temp[[i]] <- tempfile()
					}
					temp.file      <- unlist(temp.file)
					out.files.temp <- unlist(out.files.temp)
					groups <- sort(sample(1:parallel.groups, size=length(query),replace=T))
					query.groups <- lapply(X=c(1:parallel.groups),FUN=function(x){query[which(groups==x)]})
					for(i in 1:parallel.groups){
						Biostrings::writeXStringSet(x = query.groups[[i]], filepath=temp.file[[i]], append=F, format="fasta")
					}
				}
				rm(query.obj.temp)
			} else{
				rm(query.obj.temp)
				delete.query <- F
			}
		} else {
			query.path <- tempfile()
			# Increases time limit for downloading files to 1000 seconds.
			options(timeout=1000)
			conn <- utils::download.file(url=query, destfile=query.path)
			query.obj.temp        <- Biostrings::readDNAStringSet(query.path)
			names(query.obj.temp) <- mgsub(c(" ",","),c("_","_"),names(query.obj.temp))
			names(query.obj.temp) <- substring(names(query.obj.temp),first=1,last=50)
			if(is.null(parallel.groups)){
				Biostrings::writeXStringSet(x = query.obj.temp, filepath=query.path, append=F, format="fasta")
				delete.query <- T
			} else {
				temp.file <- list()
				length(temp.file) <- parallel.groups
				out.files.temp <- temp.file
				for(i in 1:parallel.groups){
					temp.file[[i]]      <- tempfile()
					out.files.temp[[i]] <- tempfile()
				}
				temp.file      <- unlist(temp.file)
				out.files.temp <- unlist(out.files.temp)
				groups <- sort(sample(1:parallel.groups, size=length(query),replace=T))
				query.groups <- lapply(X=c(1:parallel.groups),FUN=function(x){query[which(groups==x)]})
				for(i in 1:parallel.groups){
					Biostrings::writeXStringSet(x = query.groups[[i]], filepath=temp.file[[i]], append=F, format="fasta")
				}
			}
			rm(query)
		}
	}
	#### Create a path to hold a temporary output table
	if(is.null(table.out)){
		output.path  <- tempfile()
		delete.table <- T
	} else {
		output.path  <- table.out
		delete.table <- F
	}
	#### Vector of file extensions that should be present if a local database has been created to query against.
	DB.extensions        <- c(".nog",".nsq",".nhr",".nin")
	### A vectory of filenames that should be present if local NCBI database exists for the subject sequences.
	expected.db.files    <- paste(subject.path,DB.extensions,sep="")
	#### Run makeBlastDB to make a blast database if any of expected.db.files do not exist
	if(!all(file.exists(expected.db.files))){
		makeBlastDB(makeblastdb.exe.path,subject.path)
	}
	### If num.threads argument is set to "max", set num.threads equal to the number of cores available.
	if(num.threads=="max"){
		num.threads <- parallel::detectCores()
	}
	### Run blast!!
	if(is.null(parallel.groups)){
		command <- paste(blast.exe.path,"-db",subject.path,"-query",query.path,"-out",output.path,"-evalue",eval,"-outfmt",output.format,"-max_target_seqs",max.targets.per.query,"-max_hsps",max.matches.per.target,"-num_threads",num.threads,other.args)
		system(command,wait=T)
		# Need to find a way to check if the analysis is complete before doing things from here onward.
		result <- data.table::fread(output.path)
	} else {
		command.all <- list(); length(command.all) <- parallel.groups
		for(i in 1:parallel.groups){
			command.all[[i]] <- paste(blast.exe.path,"-db",subject.path,"-query",temp.file[[i]],"-out",out.files.temp[[i]],"-evalue",eval,"-outfmt",output.format,"-max_target_seqs",max.targets.per.query,"-max_hsps",max.matches.per.target,"-num_threads",num.threads,other.args)
		}
		### Run all commands in command.all, which are separated by " & " to allow for parallelization.
		command <- gsub(" & $","",paste0(paste0(unlist(command.all)," & "),collapse=""))
		system(command, wait=T)
		### Read each of the output hit tables and then merge them into one.
		result <- do.call(rbind,lapply(out.files.temp,FUN=data.table::fread))
	}
	if(output.format==6){
		colnames(result) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
	}
	if(!is.null(table.out)){
		all.matches <- write.table(x=result,file=table.out,sep="\t",quote=F,col.names=T,row.names=F)
	}
	### Delete the temporary files
	if(delete.subject){
		file.remove(subject.path)
		file.remove(expected.db.files)
	}
	if(delete.query){
		file.remove(query.path)
	}
	if(delete.table){
		file.remove(output.path)
	}
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
#' Thamnophis.sirtalis_exome <- get.seqs.from.gff(input.seqs=Thamnophis.sirtalis_genome.path,input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp)
#' 
#' ### Use tblastx to return up to 50 matches for each of exon (only first two exons in this example) of Thamnophis exome in Crotalus horridus genome
#' test.query    <- Thamnophis.sirtalis_exome[1:2]
#' test.subject  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/625/485/GCA_001625485.1_ASM162548v1/GCA_001625485.1_ASM162548v1_genomic.fna.gz"
#' test.50hits   <- blast(method="tblastx",subject=test.subject,query=test.query)

#' @title Make Local Blast Database
#' 
#' Wrapper for the NCBI makeblastdb command. This function makes an NCBI blast database (which is needed to run blast against).
#' 
#' @param makeblastdb.path Path to makeblastdb (included in the bin folder of blast)
#' @param subject.path Path to the file containing the sequences to include in the database
#' @return Creates a blast database file.
#' @export makeBlastDB
makeBlastDB  <- function(makeblastdb.path,subject.path){
	### pastes the parts into a character string that can be executed in terminal
	command  <- paste(makeblastdb.path,"-in",subject.path,"-parse_seqids -dbtype nucl")
	### calls terminal to execute the character string "command"
	system(command)
}

#' @title Check Executable
#' 
#' Checks if the input path is executable.
#' 
#' @param exe.path Character string with input path to check.
#' @return Integer 0 if exe.path is executable, otherwise 1.
#' @export check.if.executable
check.if.executable <- function(exe.path){
	system(paste("command -v",paste0("'",exe.path,"'"),">/dev/null 2>&1 || { echo >&2 ''; exit 1; }"))
}

