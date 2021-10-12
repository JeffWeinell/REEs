.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(data.table)
library(Biostrings)
library(igraph)
library(REEs)

args = commandArgs(trailingOnly=TRUE)

#workdir     <- args[1]
allseqs.path <- args[1]
matches.path <- args[2]
targetNames  <- args[3] # text string that is present in the names of all of the target sequences; = "WeinellEntry" for SnakeCap sequences
outDir       <- paste0(tools::file_path_sans_ext(matches.path),"_out")
# create outDir if it doesnt already exists
system(sprintf("[ ! -d '%s' ] && mkdir '%s' ",outDir,outDir))

# Read in query=subject sequences
seqs    <- Biostrings::readDNAStringSet(allseqs.path)
targets <- seqs[grep(targetNames,names(seqs))]

# read in blast hit table
match.data <- data.table::fread(matches.path, sep = "\t", header = F, stringsAsFactors = FALSE)
headers    <- c("qName", "tName", "pident", "matches", "misMatches", "gapopen","qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")
data.table::setnames(match.data, headers)
# filter self hits in which the query contig = subject contig
match.data <- match.data[match.data$qName!=match.data$tName,]
edgeMat    <- as.matrix(match.data[,c("qName","tName")])
hits.graph <- igraph::graph_from_edgelist(edgeMat,directed=F)
# names of all sequences that occur somewhere in the hit table
uniqueContigNamesHits <- unique(c(unique(edgeMat[,"qName"]), unique(edgeMat[,"tName"])))
# names of all sequences were searched
uniqueContigNamesAll  <- unique(names(seqs))
### Find connected components of hits.graph. Returns a list with length 3; the first entry is a vector indicating which group each sequence belongs to, the second indicates how many sequences are in each group, and the third indicates how many groups exist.
hits.clusters <- igraph::components(hits.graph)
# An array with length equal to the number of groups; the ith entry is a vector with the names of sequences in the ith group
hits.groups             <- igraph::groups(hits.clusters) # Each group is expected (in theory) to contain contigs of homologous loci
seqs.GroupMembership.df <- data.frame(SequenceName=names(hits.clusters$membership), Group=unname(hits.clusters$membership))

#### save some tables
edgeMat.path   <- file.path(outDir,paste0(basename(tools::file_path_sans_ext(matches.path)),".EdgeMatrix"))
GroupMemb.path <- file.path(outDir,paste0(basename(tools::file_path_sans_ext(matches.path)),".BlastHitGroups"))
write.table(edgeMat,edgeMat.path,col.names=F,row.names=F,sep="\t")
write.table(seqs.GroupMembership.df,GroupMemb.path,col.names=F,row.names=F,sep="\t")

#### Sort groups into three sets: (1) groups with a single reference target; (2) groups with multiple reference targets; (3) groups without a reference target, which correspond to bycatch loci.
# Groups that contain at least one target sequence
TargetGroupsIDs       <- unique(seqs.GroupMembership.df$Group[grep(targetNames,seqs.GroupMembership.df$SequenceName)])
# Number of targets in each group containing at least one target
TargetsPerTargetGroup <- table(seqs.GroupMembership.df$Group[grep(targetNames,seqs.GroupMembership.df$SequenceName)])
# Single-target groups; a list of character vectors, each with the names of sequences belonging to a group. Every group includes one target sequence.
SingleTargetGroups    <- hits.groups[names(TargetsPerTargetGroup)[which(TargetsPerTargetGroup==1)]]
# Multi-target groups; a list of character vectors, each with the names of sequences belonging to a group. Every group includes at least two different target sequences.
MultiTargetGroups    <- hits.groups[names(TargetsPerTargetGroup)[which(TargetsPerTargetGroup>1)]]
# Bycatch groups; a list of character vectors, each with the names of sequences belonging to a group. Group do not contain a target sequence.
BycatchGroups <- hits.groups[-sort(TargetGroupsIDs)] # Most of the Bycatch groups only contain two sequences.

#### Create directories for fasta files of single-target, multi-target, and bycatch groups
singleDir  <- file.path(outDir,"SingleTarget_unaligned")
multiDir   <- file.path(outDir,"MultiTarget_unaligned")
bycatchDir <- file.path(outDir,"Bycatch_unaligned")
system(sprintf("[ ! -d '%s' ] && mkdir '%s' ",singleDir,singleDir))
system(sprintf("[ ! -d '%s' ] && mkdir '%s' ",multiDir,multiDir))
system(sprintf("[ ! -d '%s' ] && mkdir '%s' ",bycatchDir,bycatchDir))

#### Write fasta files for each group to one of three directories (SingleTargets, MultiTargets, Bycatch)
# Single Target Groups
for(i in 1:length(SingleTargetGroups)){
	dna.group.i  <- seqs[SingleTargetGroups[[i]]]
	targetName.i <- names(dna.group.i)[grep(targetNames,names(dna.group.i))]
	writeXStringSet(dna.group.i,file.path(singleDir,paste0(targetName.i,".fa")))
}

# Multi Target Groups
for(i in 1:length(MultiTargetGroups)){
	dna.group.i   <- seqs[MultiTargetGroups[[i]]]
	#targetNames.i <- paste(names(dna.group.i)[grep(targetNames,names(dna.group.i))],collapse="_")
	locusName.i <- paste0("MultiTargetLocus",i)
	writeXStringSet(dna.group.i,file.path(multiDir,paste0(locusName.i,".fa")))
}

# Bycatch Groups
minGroupSize=20
BycatchGroups.use <- BycatchGroups[lengths(BycatchGroups)>=minGroupSize]
for(i in 1:length(BycatchGroups.use)){
	dna.group.i <- seqs[BycatchGroups.use[[i]]]
	locusName.i <- paste0("BycatchLocus",i)
	writeXStringSet(dna.group.i,file.path(bycatchDir,paste0(locusName.i,".fa")))
}

############
if(FALSE){
	# Moderate or strong matches
	filt.data <- match.data[match.data$matches > 40 & match.data$evalue <= 0.05 & match.data$bitscore >= 50,]
	# Remove weak matches
	matches.removed <- match.data[match.data$matches <= 40 | match.data$evalue > 0.05 | match.data$bitscore < 50,]

	# New column that holds the contig name and target name for each match
	# QuerySubject <- paste(filt.data$qName,filt.data$tName)
	filt.data[,QuerySubject:=paste(filt.data$qName,filt.data$tName)]
	filt.df      <- filt.data[match(unique(filt.data$QuerySubject),filt.data$QuerySubject),]

	# Filtered hit table for the set of contigs that match a single target locus. If the contigs are large they may contain multiple targets though.
	targets.filt.df <- filt.df[grep(targetNames,filt.df$qName),]

	filt.df.SingleMatches <- filt.df[match(names(which(table(filt.df$qName)==1)),filt.df$qName),]
	# Filtered hit table for the set of contigs that match multiple target loci. Not sure what to do with these yet
	filt.df.MultipleMatches <- filt.df[filt.df$qName %in% names(which(table(filt.df$qName)>1)),]
}

if(FALSE){
	workdir <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/"
	procdir <- file.path(workdir,"/Processed_Samples/")
	setwd(procdir)
	# paths to consensus contig fasta files (consensus among haplocontigs of an individual)
	consensusContigs.paths <- list.files(pattern='.+_consensus-contigs-dipspades.fa',recursive=T,full.names=T) # Change pattern to ".+_genomic.fna$" when dealing with genomes.
	targets.path           <- file.path(workdir,'Weinell_TargetLoci_Snakes_Final_18April2019.txt')             # Change to file.path(workdir,'OrthologAlignments_Consensus.fa') for genomes
	# paths to hit tables
	matches.paths <- list.files(pattern='.+_match.txt$',recursive=T, full.names=T)
	
	for(i in 1:length(matches.paths)){
		SampleName            <- basename(dirname(matches.paths[i]))
		match.data <- data.table::fread(matches.paths[i], sep = "\t", header = F, stringsAsFactors = FALSE)
		headers    <- c("qName", "tName", "pident", "matches", "misMatches", "gapopen","qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")
		data.table::setnames(match.data, headers)
		# Moderate or strong matches
		filt.data <- match.data[match.data$matches > 40 & match.data$evalue <= 0.05 & match.data$bitscore >= 50,]
		# Remove weak matches
		matches.removed <- match.data[match.data$matches <= 40 | match.data$evalue > 0.05 | match.data$bitscore < 50,]
		
		# Read in assembled contigs
		contigs <- Biostrings::readDNAStringSet(consensusContigs.paths[i])
		names(contigs) <- gsub(" .+","",names(contigs))
		# Read in target loci
		query.seqs <- Biostrings::readDNAStringSet(targets.path)
		
		# Adds a column that will hold direction information
		filt.data[,qDir:= "+"]
		filt.data$qDir[filt.data$tStart > filt.data$tEnd] <- "-"
		
		# New column that holds the contig name and target name for each match
		# ContigTarget <- paste(filt.data$qName,filt.data$tName)
		filt.data[,ContigTarget:=paste(filt.data$qName,filt.data$tName)]
		filt.df <- filt.data[match(unique(filt.data$ContigTarget),filt.data$ContigTarget),]
		
		# Filtered hit table for the set of contigs that match a single target locus. If the contigs are large they may contain multiple targets though.
		filt.df.SingleMatches <- filt.df[match(names(which(table(filt.df$qName)==1)),filt.df$qName),]
		# Filtered hit table for the set of contigs that match multiple target loci. Not sure what to do with these yet
		filt.df.MultipleMatches <- filt.df[filt.df$qName %in% names(which(table(filt.df$qName)>1)),]
		
		# Contigs that strongly match only one target
		contigs.single        <- contigs[pmatch(filt.df.SingleMatches$qName,names(contigs))]
		newNames.single       <- paste0(names(contigs.single),"_|_",filt.df.SingleMatches$tName,"_|_",SampleName)
		names(contigs.single) <- newNames.single
		# Write the set of contigs that match a single target locus, and include the sample name, contig name, and target name in the sequence name
		Biostrings::writeXStringSet(contigs.single,file.path(dirname(matches.paths[i]),paste0(SampleName,"_contigs-SingleTargetMatches.fa")))
		# Write the filtered Hit table for the contigs that match a single target
		write.table(filt.df.SingleMatches,file.path(dirname(matches.paths[i]),paste0(SampleName,"_Hits-SingleTargetMatches.tsv")),quote=F,sep="\t",col.names=T,row.names=F)
		
		# Contigs that moderately or strongly match multiple target
		contigs.multi        <- contigs[match(filt.df.MultipleMatches$qName,names(contigs))]
		newNames.multi       <- paste0(names(contigs.multi),"_|_",SampleName)
		names(contigs.multi) <- newNames.multi
		# Write the set of contigs that match multiple target loci, and include the sample name and contig name in the sequence name
		Biostrings::writeXStringSet(contigs.multi,file.path(dirname(matches.paths[i]),paste0(SampleName,"_contigs-MultipleTargetMatches.fa")))
		# Write the filtered Hit table for the contigs that match multiple targets
		write.table(filt.df.MultipleMatches,file.path(dirname(matches.paths[i]),paste0(SampleName,"_Hits-MultipleTargetMatches.tsv")),quote=F,sep="\t",col.names=T,row.names=F)
		
		# Contigs that only weakly match one or more targets
		if(any(!matches.removed$qName %in% filt.data$qName)){
			contigs.weakMatches <- contigs[match(unique(matches.removed$qName[which(!matches.removed$qName %in% filt.data$qName)]),names(contigs))]
			newNames.weak       <- paste0(names(contigs.weakMatches),"_|_",SampleName)
			# Write the set of contigs that had only weak matche(s) to targets, and include the sample name and contig name in the sequence name
			Biostrings::writeXStringSet(contigs.weakMatches,file.path(dirname(matches.paths[i]),paste0(SampleName,"_contigs-WeakTargetMatchesOnly.fa")))
			# write the hit table containing only weak matches
			write.table(matches.removed,file.path(dirname(matches.paths[i]),paste0(SampleName,"_Hits-WeakMatches.tsv")),quote=F,sep="\t",col.names=T,row.names=F)
		}
		
		# Contigs that do not match any targets
		if(any(!names(contigs) %in% match.data$qName)){
			contigs.unmatched  <- contigs[!names(contigs) %in% match.data$qName]
			newNames.unmatched <- paste0(names(contigs.unmatched),"_|_",SampleName)
			# Write the set of contigs that did not match any of the targets
			Biostrings::writeXStringSet(contigs.unmatched,file.path(dirname(matches.paths[i]),paste0(SampleName,"_contigs-unmatched.fa")))
		}
	
	}
	#####
	## Is it possible to write the code above in bash/linux, and then merge with targetMatching.sh script?
	#####
}


