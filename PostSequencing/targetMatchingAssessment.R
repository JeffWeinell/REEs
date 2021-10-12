.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(data.table)
library(Biostrings)

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
	ContigTarget <- paste(filt.data$qName,filt.data$tName)
	filt.data[,ContigTarget:=paste(filt.data$qName,filt.data$tName)]
	filt.df      <- filt.data[match(unique(filt.data$ContigTarget),filt.data$ContigTarget),]
	
	# Filtered hit table for the set of contigs that match a single target locus. If the contigs are large they may contain multiple targets though.
	filt.df.SingleMatches   <- filt.df[match(names(which(table(filt.df$qName)==1)),filt.df$qName),]
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




