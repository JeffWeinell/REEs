.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(REEs)

# workdir <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/"
# targets <- file.path(workdir,'Weinell_TargetLoci_Snakes_Final_18April2019.txt')
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
	stop("at least one argument must be supplied to indicate the path to the target sequences")
} 
# else if(length(args)==1){
# 	args[2] <- getwd()
# 	args[3] <- 1
# 	args[4] <- NULL
# } else if(length(args)==2){
# 	args[3] <- 1
# 	args[4] <- NULL
# } else if(length(args)==3){
# 	args[4] <- NULL
# }

res     <- lociAlignment04(targets.path=args[1],workdir=args[2],istart=args[3],iend=args[4])

if(FALSE){
	library(data.table)
	library(Biostrings)
	library(REEs)
	workdir <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/"
	procdir <- file.path(workdir,"/Processed_Samples")
	setwd(procdir)
	
	# paths to consensus sequence fasta files
	consensusContigs.paths <- list.files(pattern='.+_consensus-contigs-dipspades.fa',recursive=T,full.names=T)
	targets.path <- file.path(workdir,'Weinell_TargetLoci_Snakes_Final_18April2019.txt')
	# paths to hit tables
	# matches.paths <- list.files(pattern='.+_match.txt$',recursive=T, full.names=T)
	
	# Read in assembled contigs
	# contigs <- Biostrings::readDNAStringSet(consensusContigs.paths)
	# Read in target loci
	query.seqs <- Biostrings::readDNAStringSet(targets.path)
	
	#####
	# paths to consensus contigs that strongly match only one target
	SingleTargetContigs.paths <- list.files(pattern='.+_contigs-SingleTargetMatches.fa',recursive=T,full.names=T)
	# DNAStringSet holding the consensus contigs that strongly match only one target, for all of the species sampled
	SingleTargetContigs <- Biostrings::readDNAStringSet(SingleTargetContigs.paths)
	dir.check.create(file.path(workdir,"/SingleTargetContigs"))
	setwd(file.path(workdir,"/SingleTargetContigs"))
	AlignmentDir1 <- file.path(workdir,"/SingleTargetContigs/UntrimmedAlignments_SamplesOnly/")
	AlignmentDir2 <- file.path(workdir,"/SingleTargetContigs/UntrimmedAlignments_SamplesAndTarget/")
	ConsensusDir  <- file.path(workdir,"/SingleTargetContigs/AlignmentConsensus_SamplesAndTarget/")
	dir.check.create(AlignmentDir1)
	dir.check.create(AlignmentDir2)
	dir.check.create(ConsensusDir)
	unresolvedDir <- file.path(workdir,"/SingleTargetContigs/UnresolvedHomologs_unaligned/")
	paralogsDir   <- file.path(workdir,"/SingleTargetContigs/PutativeParalogs_unaligned/")
	orthologsDir  <- file.path(workdir,"/SingleTargetContigs/PutativeOrthologs_unaligned/")
	dir.check.create(unresolvedDir)
	dir.check.create(paralogsDir)
	dir.check.create(orthologsDir)
	# Directory where sequences will be written if the target was only sequenced in one individual.
	underSampledDir <- file.path(workdir,"/SingleTargetContigs/undersampledTargets_unaligned/")
	dir.check.create(underSampledDir)
	
	# For each target, write a separate fasta file containing the set of species' contigs that match it
	contigInfo.df <- as.data.frame(do.call(rbind,strsplit(names(SingleTargetContigs),split="_|_",fixed=T)))
	colnames(contigInfo.df) <- c("contig","target","sample")
	uniqueTargets <- unique(contigInfo.df$target)
	# For each target and each individual, the number of sequence contigs
	SampleContigsPerTarget <- table(contigInfo.df[,c("target","sample")])
	IndvsPerTarget <- SampleContigsPerTarget
	IndvsPerTarget[SampleContigsPerTarget!=0] <- 1
	
	print(sprintf("%s target loci with %s–%s contigs/target and %s–%s individuals/target",length(uniqueTargets), min(rowSums(SampleContigsPerTarget)),max(rowSums(SampleContigsPerTarget)),min(rowSums(IndvsPerTarget)),max(rowSums(IndvsPerTarget))))
	for(i in 1:length(uniqueTargets)){
	#for(i in 12:13){
		# DNAStringSet of sample contigs for the ith target
		contigs.i      <- SingleTargetContigs[contigInfo.df$target==uniqueTargets[i]]
		# Number of times that each individual occurs in the sequence set. Would only be once if only orthologs present (assuming consensus haplocontigs are used).
		countSamples   <- table(gsub(".+_\\|_.*_\\|_","",names(contigs.i)))
		# Jump to next target locus if only one individual sequenced
		if(length(countSamples)==1){
			Biostrings::writeXStringSet(contigs.i,file.path(underSampledDir,paste0("/",uniqueTargets[i],"_undersampled.fa")))
			next
		}
		# Jump to next target locus if more than 100 possible alternative alignments
		if(prod(countSamples)>100){
			Biostrings::writeXStringSet(contigs.i,file.path(unresolvedDir,paste0("/",uniqueTargets[i],"_HomologyUnresolved.fa")))
			next
		}
		# DNAStringSet of sample contigs and the target for the ith target
		seqs.i         <- c(query.seqs[uniqueTargets[i]],contigs.i)
		# A list of vectors. Each vector holds contig names for putative homologs of the target sequence, one of which is assumed to be an ortholog and the others paralogs to the target sequence.
		alternativeContigs.list <- lapply(names(countSamples),function(x) names(contigs.i)[grep(x,names(contigs.i))])
		# Data frame with each row holding the contig names to include in an alignment. Each row has a different set of contig names, and each column represents an individual. The best-scoring alignment will be saved and assumed to contain orthologous contigs.
		alternativeContigs.df <- expand.grid(alternativeContigs.list)
		alignments.list <- list(); length(alignments.list) <- nrow(alternativeContigs.df)
		pdist.list      <- list(); length(pdist.list)      <- nrow(alternativeContigs.df)
		for(j in 1:nrow(alternativeContigs.df)){
		#for(j in 1:4){
			contigs.ij <- seqs.i[as.character(unlist(alternativeContigs.df[j,]))]
			# including the target
			seqs.ij    <- c(query.seqs[uniqueTargets[i]],contigs.ij)
			# alignment for the jth alternative set of contigs for the ith target locus
			aligned.ij <- Biostrings::DNAStringSet(REEs::mafft(seqs.ij, param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6"))
			# mean pairwise distance of aligned.ij; for now this will be a measure of the 
			meanpdist.ij <- mean(ape::dist.dna(ape::as.DNAbin(aligned.ij),pairwise.deletion=T),na.rm=T)
			alignments.list[[j]] <- aligned.ij
			pdist.list[[j]] <- meanpdist.ij
			print(sprintf("alternative contigs alignment %s of %s determined for target '%s' (%s of %s)",j,nrow(alternativeContigs.df),uniqueTargets[i],i,length(uniqueTargets)))
		}
		# Names of the putatively orthologous contigs
		bestContigs.names   <- as.character(unlist(alternativeContigs.df[which(unlist(pdist.list)==min(unlist(pdist.list)))[1],]))
		# The alignment that includes the best contigs and the reference target
		alignedContigsTarget.best <- alignments.list[[which(unlist(pdist.list)==min(unlist(pdist.list)))[1]]]
		# Aligns the best contigs without including the target reference
		contigs.best        <- contigs.i[bestContigs.names]
		alignedContigs.best <- Biostrings::DNAStringSet(REEs::mafft(contigs.best, param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6"))
		# Writing the best alignment (with and without the reference target) for the ith target
		Biostrings::writeXStringSet(alignedContigs.best,file.path(AlignmentDir1,paste0("/",uniqueTargets[i],".fa")))
		Biostrings::writeXStringSet(alignedContigsTarget.best,file.path(AlignmentDir2,paste0("/",uniqueTargets[i],".fa")))
		# Determining and then writing the consensus sequence of the best alignment that includes contigs and the reference target
		consensus.of.alignedContigsTarget.best <- Biostrings::DNAStringSet(Biostrings::consensusString(Biostrings::DNAStringSet(gsub("-","N",alignedContigsTarget.best)),ambiguityMap=IUPAC_CODE_MAP))
		names(consensus.of.alignedContigsTarget.best) <- paste0(uniqueTargets[i],"_consensus")
		Biostrings::writeXStringSet(consensus.of.alignedContigsTarget.best,file.path(ConsensusDir,paste0("/",uniqueTargets[i],"_consensus.fa")))
		# Writing unaligned putative orthologs of the ith target to a fasta file
		Biostrings::writeXStringSet(contigs.best,file.path(orthologsDir,paste0("/",uniqueTargets[i],"_orthologs.fa")))
		# Writing unaligned putative paralogs of the ith target to a fasta file
		if(nrow(alternativeContigs.df) > 1){
			paralogs.i  <- contigs.i[setdiff(names(contigs.i),bestContigs.names)]
			paralogInfo <- as.data.frame(do.call(rbind,strsplit(names(paralogs.i),split="_|_",fixed=T)))
			paralogNames <- sapply(1:nrow(paralogInfo),function(x) paste(data.frame(V1=paralogInfo[x,1],V2=paste0(paralogInfo[x,2],"-putative-paralog"),V3=paralogInfo[x,3]),collapse="_|_"))
			names(paralogs.i) <- paralogNames
			Biostrings::writeXStringSet(paralogs.i,file.path(paralogsDir,paste0("/",uniqueTargets[i],"_paralogs.fa")))
		}
		print(sprintf("ortholog alignment of %s (target %s of %s) complete!",uniqueTargets[i],i,length(uniqueTargets)))
	}
}

