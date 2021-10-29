.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(data.table)
library(igraph)

args       <- commandArgs(trailingOnly=TRUE)
workdir    <- args[1]
outpath    <- args[2]
i.start    <- as.numeric(args[3])
i.end      <- as.numeric(args[4])

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
write.table(groups.df, outpath, col.names=T, row.names=F, sep="\t")


