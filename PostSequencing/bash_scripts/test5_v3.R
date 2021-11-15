.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(data.table)
library(gtools)

args       <- commandArgs(trailingOnly=TRUE)
workdir    <- args[1] # workdir="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/"
outpath    <- args[2]
i.start    <- as.numeric(args[3])
i.end      <- as.numeric(args[4])
mode       <- args[5]

if(mode=="B"){
	library(igraph)
	groupsDir  <- file.path(workdir,"groups")
	sampleskey <- file.path(workdir,"allsamples_key.txt")
	
	print(sprintf("Preparing to construct group matrix for group matrices in '%s'", groupsDir))
	
	# filepaths to all of the pairwise group matrices
	filepaths0 <- sort(list.files(groupsDir, full.names=T))
	
	sampleskey.df     <- read.table(sampleskey,header=F,sep=" ")
	samples.SeqCap    <- which(sampleskey.df[,3]=="SeqCap")
	samples.genome    <- which(sampleskey.df[,3]=="genome")
	samples.targetRef <- which(sampleskey.df[,3]=="targetRef")
	
	# samples 1–35 = sequence capture samples
	# samples 35–65 = genomes
	# sample 66 = reference targets
	
	filepaths.capVScap        <- filepaths0[-sort(unique(unlist(lapply(paste0("sample",c(samples.genome,samples.targetRef),"_"), function(x){grep(x,filepaths0)}))))]
	filepaths.capVSgenome     <- filepaths0[sort(intersect(unlist(lapply(paste0("sample",samples.SeqCap,"_"), function(x){grep(x,filepaths0)})), unlist(lapply(paste0("sample",samples.genome,"_"), function(x){grep(x,filepaths0)}))))]
	filepaths.capVStargets    <- filepaths0[sort(intersect(unlist(lapply(paste0("sample",samples.SeqCap,"_"), function(x){grep(x,filepaths0)})), unlist(lapply(paste0("sample",samples.targetRef,"_"), function(x){grep(x,filepaths0)}))))]
	filepaths.genomeVSgenome  <- filepaths0[-sort(unique(unlist(lapply(paste0("sample",c(samples.SeqCap,samples.targetRef),"_"), function(x){grep(x,filepaths0)}))))]
	filepaths.targetsVSgenome <- filepaths0[sort(intersect(unlist(lapply(paste0("sample",samples.genome,"_"), function(x){grep(x,filepaths0)})), unlist(lapply(paste0("sample",samples.targetRef,"_"), function(x){grep(x,filepaths0)}))))]
	
	# filepaths to the pairwise group matrices that should be used for forming putative orthology groups
	filepaths  <- c(filepaths.capVScap, filepaths.capVStargets)
	# filepaths  <- c(filepaths.capVScap, filepaths.capVStargets, filepaths.capVSgenome, filepaths.targetsVSgenome)
	
	groups.df <- NULL
	for(i in i.start:i.end){
		groups.df.i     <- data.table::fread(filepaths[i],header=F, sep="\t")
		groups.df.i[,2] <- paste(i,unlist(groups.df.i[,"V2"]), sep="_")
		groups.df       <- rbind(groups.df,groups.df.i)
		groups.g        <- igraph::graph_from_data_frame(groups.df)
		Membership      <- igraph::clusters(groups.g)$membership
		samplerows      <- grep("sample",names(Membership))
		groups.df       <- data.table(V1=names(Membership)[samplerows], V2=unname(Membership)[samplerows])
		print(i)
	}
	data.table::setnames(groups.df,c("sequence","group"))
	summary.df   <- data.frame(numFiles=i, rows.groups.df=dim(groups.df)[1], NumGroups.groups.df=length(unique(unlist(groups.df[,2]))))
	print(summary.df)
	write.table(groups.df, outpath, col.names=T, row.names=F, quote=F, sep="\t")
}

# Concatenates the set of partially combined group matrices into a fully combined group matrix
if(mode=="C"){
	# Code not yet implemented
}

if(mode=="D"){
	#allseqsMat.path <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/allseqs_key2.txt"
	#groupsMat       <- data.table::fread("/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupsMatrix_capVScap_capVStargets_v3.txt",sep="\t",header=F)
	groupsMat       <- args[1]
	outpath         <- args[2]
	allseqsMat.path <- args[4]
	# Need to capture allseqsMat.path and groupsMat from the SH file!
	# empty lists to fill
	groupName    <- unique(unname(unlist(groupsMat[,"V1"])))
	nSamples     <- list(); length(nSamples) <- length(groupName)
	nSeqs        <- list(); length(nSeqs) <- length(groupName)
	nTargs       <- list(); length(nTargs) <- length(groupName)
	targSeqNames <- list(); length(targSeqNames) <- length(groupName)
	targNames    <- list(); length(targNames) <- length(groupName)
	
	for(i in 1:length(groupName)){
		name.temp     <- groupName[i]
		mat.temp      <- as.matrix(groupsMat[which(unname(unlist(groupsMat[,"V1"]))==name.temp),])
		seqs.temp     <- mat.temp[,2]
		nSeqs[[i]]    <- nrow(mat.temp)
		samples.temp  <- unique(gsub("_.+","",seqs.temp))
		nSamples[[i]] <- length(samples.temp)
		if(any(samples.temp=="sample66")){
			targSeqNames.temp <- seqs.temp[grep("sample66_",seqs.temp)]
			targSeqNames[[i]] <- paste(targSeqNames.temp,collapse="|")
			targNames.temp    <- system(sprintf("awk -v x='%s' '{if($1 ~ x) {print $2}}' '%s'", targSeqNames[[i]], allseqsMat.path),intern=T)
			targNames[[i]]    <- paste(targNames.temp,collapse="|")
			nTargs[[i]]       <- length(targSeqNames.temp)
		} else {
			targSeqNames[[i]] <- "NA"
			targNames[[i]]    <- "NA"
			nTargs[[i]]       <- 0
		}
		#nSamples     <- unlist(nSamples)
		#nSeqs        <- unlist(nSeqs)
		#nTargs       <- unlist(nTargs)
		#targSeqNames <- unlist(targSeqNames)
		#targNames    <- unlist(targNames)
		result       <- data.frame(group=groupName[i],N_UniqueSamples=nSamples[[i]],N_Seqs=nSeqs[[i]],N_targetRef=nTargs[[i]],targetRef_seqNames=targSeqNames[[i]],targetRef_Names=targNames[[i]])
		#outpath     <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupseqs/N_UniqueSamples_Seqs_TargRefs_PerGroup_WithHeader_v2.txt"
		if(i==1){
			APPEND=F
			HEADER=T
		} else {
			APPEND=T
			HEADER=F
		}
		write.table(result,file=outpath,col.names=HEADER,row.names=F,sep="\t",quote=F,append=APPEND)
	}
}








