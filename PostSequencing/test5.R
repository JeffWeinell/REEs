.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(data.table)
library(igraph)

args       <- commandArgs(trailingOnly=TRUE)
groupsDir  <- args[1]
outpath    <- args[2]

print(sprintf("Preparing to construct group matrix for group matrices in '%s'", groupsDir))

filepaths <- sort(list.files(groupsDir, full.names=T))
groups.df <- NULL
for(i in 1:length(filepaths)){
	groups.df.i     <- data.table::fread(filepaths[i],header=F, sep="\t")
	groups.df.i[,2] <- paste(i,unlist(groups.df.i[,"V2"]), sep="_")
	groups.df       <- rbind(groups.df,groups.df.i)
	print(i)
}
data.table::setnames(groups.df,c("sequence","group"))

groups.g      <- igraph::graph_from_data_frame(groups.df)
Membership    <- igraph::clusters(groups.g)$membership
Membership.df <- data.frame(sequence=names(Membership), group=unname(Membership))
summary.df    <- data.frame(numFiles=i, rows.groups.df=dim(groups.df)[1], NumGroups.groups.df=length(unique(unlist(groups.df[,2]))), rows.Membership.df=dim(groups.df)[1], NumGroups.Membership.df=length(unique(unlist(Membership.df[,2]))))
summary.df
write.table(Membership.df, outpath, col.names=T, row.names=F, sep="\t")
