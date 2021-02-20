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
#' XStringSet object (DNAStringSet or AAStringSet). If a XStringSet object is used, then a temporary fasta file is created to hold the sequences to search within.
#' If DNAStringSet or URL character string, a local copy will be saved as a temporary file.
#' @param query Path to the fasta file containing the query sequences to find matches for. If a XStringSet object is used, then a temporary fasta file is created to hold the query sequences.
#' @param table.out Path where to save the output table (including filename).
#' @param eval Number to use for the expect value (e-value); default is 1e-5.
#' @param output.format Integer indicating which format to use for the output table of matches (default 6).
#' @param max.targets.per.query Integer indicating the maximum number of targets (subject contigs) containing a match per query sequence. This sets the BLAST's max_target_seqs argument. Importantly, this argument does not control the total number of matches per query, because multiple matches may occer on the same subject contig.
#' @param max.matches.per.target Integer indicating the maximum number of matches per target (subject contig) per query. This sets BLAST's max_hsps argument. To determine the maximum total number of matches per query, multiply max.targets.per.query by max.matches.per.target.
#' @param parallel.groups Number of groups to split the query sequences into for running in parallel. Default 10. Will likely produce an error if less than the number of query sequences.
#' @param num.threads Either an integer indicating how many threads to use, or "max" (default), in which case num.threads is set to the number of cores available.
#' @param cleanup Logical indicating if temporary output files should be deleted. Default is TRUE.
#' @param other.args A character string of the form "-argument1 value1 -argument2 value2" with additional arguments and values to pass to BLAST. See BLAST manual for definitions of available arguments. Default is NULL.
#' @return Table of matches.
#' @export blast
blast <- function(blast.path="auto",method,subject,query,table.out,eval=1e-5,output.format=6,max.targets.per.query=10,max.matches.per.target=10,parallel.groups=10,num.threads="max",cleanup=T,other.args=NULL){
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
	# Test that table.out does not already exist as a file or directory.
	if(dirname(table.out) == "."){
		table.out <- paste0(getwd(),"/",table.out)
	}
	if(file.exists(table.out) | dir.exists(table.out) | !dir.exists(dirname(table.out))){
		stop("invalid table.out")
	}
	# Where subject and query files data should be held.
	output.path       <- table.out
	# Create a directory to hold output files other than output.table; this directory will be deleted upon close if cleanup = TRUE
	temp.path         <- paste0(dirname(output.path),"/",basename(tempfile(pattern="tempoutdir")))
	dircon            <- dir.check.create(temp.path)
	subject.path.temp <- paste0(temp.path,"/",basename(tempfile()))
	query.path.temp   <- paste0(temp.path,"/",basename(tempfile()))
	if(!is.null(parallel.groups)){
		out.files.temp   <- paste0(temp.path,"/",basename(paste0(tools::file_path_sans_ext(table.out),"_",c(1:parallel.groups),".tsv")))
		query.paths.all  <- paste0(query.path.temp,"_",c(1:parallel.groups))
	}
	if(is(subject,"XStringSet")){
		names(subject) <- gsub(" .+","",names(subject))
		print("preparing subject sequences")
		scon1 <- Biostrings::writeXStringSet(x = subject, filepath=subject.path.temp, append=F, format="fasta")
	} else {
		if(file.exists(subject)){
			subject.path     <- subject
			if(method %in% c("blastn","tblastn","tblastx")){
				# Check on what these do:"rpsblast", "rpstblastn"
				subject.obj.temp <- Biostrings::readDNAStringSet(subject.path)
			} else {
				if(method %in% c("blastp","deltablast","psiblast","blastx")) {
					subject.obj.temp <- Biostrings::readAAStringSet(subject.path)
				} else {
					stop("'method' argument must be 'blastn','tblastn','tblastx', 'blastp', or 'blastx'")
				}
			}
			names(subject.obj.temp) <- gsub(" .+","",names(subject.obj.temp))
			print("preparing subject sequences")
			con3         <- Biostrings::writeXStringSet(x = subject.obj.temp, filepath=subject.path.temp, append=F, format="fasta")
		} else {
			# sets time limit for downloading files to 1000 seconds
			options(timeout=1000)
			conn <- utils::download.file(url=subject, destfile=subject.path,method="auto")
			if(method %in% c("blastn","tblastn","tblastx")){
				subject.obj.temp <- Biostrings::readDNAStringSet(subject.path)
			} else{
				if(method %in% c("blastp","deltablast","psiblast","blastx")) {
					subject.obj.temp <- Biostrings::readAAStringSet(subject.path)
				} else {
					stop("'method' argument must be 'blastn','tblastn','tblastx', 'blastp', or 'blastx'")
				}
			}
			names(subject.obj.temp) <- gsub(" .+","",names(subject.obj.temp))
			print("preparing subject sequences")
			con3  <- Biostrings::writeXStringSet(x = subject.obj.temp, filepath=subject.path.temp, append=F, format="fasta")
		}
	}
	#### Make a temporary fasta file holding the query sequences if the value of the query argument is an XStringSet or URL string.
	if(is(query,"XStringSet")){
		names(query) <- mgsub(c(" ",","),c("_","_"),names(query))
		names(query) <- substring(names(query),first=1,last=50)
		if(is.null(parallel.groups)){
			print("preparing query sequences")
			qcon4 <- Biostrings::writeXStringSet(x = query, filepath=query.path.temp, append=F, format="fasta")
		} else {
			### Ensure that number of parallel groups is not more than number of query sequences.
			if(parallel.groups>length(query)){
				parallel.groups <- length(query)
			}
			groups         <- sort(sample(1:parallel.groups, size=length(query),replace=T))
			while(!all(1:parallel.groups %in% groups)){
				groups     <- sort(sample(1:parallel.groups, size=length(query),replace=T))
			}
			query.groups   <- lapply(X=c(1:parallel.groups),FUN=function(x){query[which(groups==x)]})
			print("preparing query sequences")
			for(i in 1:parallel.groups){
				qcon.temp <- Biostrings::writeXStringSet(x = query.groups[[i]], filepath=query.paths.all[i], append=F, format="fasta")
			}
		}
	} else {
		if(file.exists(query)){
			query.path     <- query
			if(method %in% c("blastn","tblastx","blastx")){
				query.obj.temp <- Biostrings::readDNAStringSet(query.path)
			} else {
				if(method %in% c("tblastn","blastp","deltablast","psiblast")){
					query.obj.temp <- Biostrings::readAAStringSet(query.path)
				} else {
					stop("'method' argument must be 'blastn','tblastn','tblastx', 'blastp', or 'blastx'")
				}
			}
			names(query.obj.temp) <- substring(mgsub(c(" ",","),c("_","_"),names(query.obj.temp)),first=1,last=50)
			if(is.null(parallel.groups)){
				print("preparing query sequences")
				qcon5 <- Biostrings::writeXStringSet(x = query.obj.temp, filepath=query.path.temp, append=F, format="fasta")
			} else {
				### Ensure that number of parallel groups is not more than number of query sequences.
				if(parallel.groups>length(query)){
					parallel.groups <- length(query)
				}
				groups <- sort(sample(1:parallel.groups, size=length(query),replace=T))
				while(!all(1:parallel.groups %in% groups)){
					groups         <- sort(sample(1:parallel.groups, size=length(query),replace=T))
				}
				query.groups <- lapply(X=c(1:parallel.groups),FUN=function(x){query[which(groups==x)]})
				print("preparing query sequences")
				for(i in 1:parallel.groups){
					qcon.temp <- Biostrings::writeXStringSet(x = query.groups[[i]], filepath=query.paths.all[i], append=F, format="fasta")
				}
			}
		} else {
			# Increases time limit for downloading files to 5000 seconds.
			options(timeout=5000)
			conn <- utils::download.file(url=query, destfile=query.path,method="auto")
			if(method %in% c("blastn","tblastx","blastx")){
				query.obj.temp <- Biostrings::readDNAStringSet(query.path)
			}
			if(method %in% c("tblastn","blastp","deltablast","psiblast")){
				query.obj.temp <- Biostrings::readAAStringSet(query.path)
			}
			query.obj.temp        <- Biostrings::readDNAStringSet(query.path)
			names(query.obj.temp) <- mgsub(c(" ",","),c("_","_"),names(query.obj.temp))
			names(query.obj.temp) <- substring(names(query.obj.temp),first=1,last=50)
			if(is.null(parallel.groups)){
				print("preparing query sequences")
				qcon6 <- Biostrings::writeXStringSet(x = query.obj.temp, filepath=query.path.temp, append=F, format="fasta")
				delete.query <- T
			} else {
				### Ensure that number of parallel groups is not more than number of query sequences.
				if(parallel.groups>length(query)){
					parallel.groups <- length(query)
				}
				groups <- sort(sample(1:parallel.groups, size=length(query),replace=T))
				while(!all(1:parallel.groups %in% groups)){
					groups         <- sort(sample(1:parallel.groups, size=length(query),replace=T))
				}
				query.groups <- lapply(X=c(1:parallel.groups),FUN=function(x){query[which(groups==x)]})
				print("preparing query sequences")
				for(i in 1:parallel.groups){
					qcon.temp <- Biostrings::writeXStringSet(x = query.groups[[i]], filepath=query.paths.all[i], append=F, format="fasta")
				}
			}
		}
	}
	#### Vector of file extensions that should be present if a local database has been created to query against.
	DB.extensions      <- c(".nog",".nsq",".nhr",".nin")
	### A vectory of filenames that should be present if local NCBI database exists for the subject sequences.
	expected.db.files  <- paste0(subject.path.temp,DB.extensions)
	#### Run makeBlastDB to make a blast database if any of expected.db.files do not exist
	if(!all(file.exists(expected.db.files))){
		print("creating blast database")
		dbcon <- makeBlastDB(makeblastdb.exe.path,subject.path.temp)
	}
	### If num.threads argument is set to "max", set num.threads equal to the number of cores available.
	if(is.null(parallel.groups)){
		if(num.threads=="max"){
			num.threads <- parallel::detectCores()
		}
	} else {
		if(num.threads=="max"){
			num.threads <- max(floor(parallel::detectCores()/parallel.groups),1)
		} else {
			num.threads <- max(floor(num.threads/parallel.groups),1)
		}
	}
	### Run blast!!
	print("about to run blast")
	if(is.null(parallel.groups)){
		command <- paste(blast.exe.path,"-db",subject.path.temp,"-query",query.path.temp,"-out",output.path,"-evalue",eval,"-outfmt",output.format,"-max_target_seqs",max.targets.per.query,"-max_hsps",max.matches.per.target,"-num_threads",num.threads,other.args)
		res <- system(command, wait=T)
		if(res == 0) {
			result <- data.table::fread(output.path)
		} else {
			stop(paste(method,"had nonzero exit status"))
		}
	} else {
		command.all <- paste(blast.exe.path,"-db",subject.path.temp,"-query",query.paths.all,"-out",out.files.temp,"-evalue",eval,"-outfmt",output.format,"-max_target_seqs",max.targets.per.query,"-max_hsps",max.matches.per.target,"-num_threads",num.threads,other.args)
		res         <- parallel::mclapply(X=command.all,FUN=system,wait=T)
		if(any(unlist(res)!=0)){
			warn.incomplete <- warning(paste("nonzero exit status for groups",which(unlist(res)!=0)))
		}
		### Now check which output files are nonempty
		if(any(file.size(out.files.temp)!=0)){
			out.files.nonempty <- out.files.temp[which(file.size(out.files.temp)!=0)]
		} else {
			stop("all output files empty")
		}
		### Read each of the output hit tables and then merge them into one.
		result <- do.call(rbind,lapply(out.files.nonempty,FUN=data.table::fread))
	}
	if(output.format==6){
		colnames(result) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
	}
	if(!is.null(table.out)){
		all.matches <- write.table(x=result,file=table.out,sep="\t",quote=F,col.names=T,row.names=F)
	}
	### Delete the folder intermediate output files
	if(cleanup==TRUE){
		system(paste("rm -R",temp.path))
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
	command  <- paste(makeblastdb.path,"-in",subject.path,"-parse_seqids -dbtype nucl -max_file_sz 4GB")
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

