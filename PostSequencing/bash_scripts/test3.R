.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(data.table)
library(igraph)
# Avoid loading Biostrings or REEs packages in this script. They use a lot of memory.

# Description:
# Reads in each pair of BLAST hit matrices (M1 & M2) in which M1 Query = M2 Subject, and M1 Subject = M2 Query; /
# uses first to columns of M1 and M2 as edge matrices, and then merges and simplifies the pair of edge matrices into a single edge matrix for samples i and j.

args       <- commandArgs(trailingOnly=TRUE)
blastdbdir <- args[1]
i.start    <- as.numeric(args[2])
i.end      <- as.numeric(args[3])

matchesDir  <- file.path(blastdbdir,"matches")
edgeMatsDir <- file.path(blastdbdir,"edgeMats")
groupsDir   <- file.path(blastdbdir,"groups")

print(sprintf("Preparing to construct edge matrix for VAR=%s-%s",i.start,i.end))

if(!dir.exists(edgeMatsDir)){
	dir.create(edgeMatsDir)
}
if(!dir.exists(groupsDir)){
	dir.create(groupsDir)
}

inpaths          <- sort(list.files(matchesDir,full.names=T,pattern="_matches.txt"))
filenames        <- basename(inpaths)
samples          <- gsub("Ssample","",gsub("Qsample","",gsub("_matches.txt","",filenames)))
samplesMat       <- do.call(rbind,strsplit(samples,split="_"))
mode(samplesMat) <-"numeric"
paths.mat        <- data.frame(path=inpaths,Qname=samplesMat[,1],Sname=samplesMat[,2])
ijdf             <- paths.mat[(paths.mat$Qname > paths.mat$Sname),]
colnames(ijdf) <- c("path","i","j")
for(x in i.start:i.end){
	i.temp <- ijdf$i[x]
	j.temp <- ijdf$j[x]
	path1  <- paths.mat$path[paths.mat$Qname == i.temp & paths.mat$Sname == j.temp]
	path2  <- paths.mat$path[paths.mat$Qname == j.temp & paths.mat$Sname == i.temp]
	# reads first two columns of hit tables
	table1 <- data.table::fread(path1, header=F, select=c(1,2))
	table2 <- data.table::fread(path2, header=F, select=c(1,2))
	table.both <- rbind(table1,table2)
	colnames(table.both) <- c("qName", "sName")
	edges.df   <- igraph::as_data_frame(igraph::simplify(igraph::graph_from_data_frame(table.both, directed=FALSE)))
	outpath.edgeMat    <- file.path(edgeMatsDir,gsub("Ssample","sample",gsub("Qsample","sample",gsub("_matches.txt","_EdgeMatrix.txt",basename(ijdf$path[x])))))
	write.table(edges.df, file=outpath.edgeMat, quote=F, sep="\t", row.names=F, col.names=T, append=F)
	edgeMat            <- as.matrix(edges.df)
	hits.clusters      <- igraph::components(igraph::graph_from_edgelist(edgeMat, directed=F))
	GroupMembership.df <- data.frame(sequence=names(hits.clusters$membership), group=unname(hits.clusters$membership))
	outpath.groups     <- file.path(groupsDir,gsub("_EdgeMatrix","_groups",basename(outpath.edgeMat)))
	write.table(GroupMembership.df, file=outpath.groups, col.names=F,row.names=F,sep="\t",quote=F)
	print(x)
}



