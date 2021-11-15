.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(Biostrings)
library(gtools)

args     <- commandArgs(trailingOnly=TRUE)
indir    <- args[1]
outdir   <- args[2]
varj     <- as.numeric(args[3])
varj2    <- as.numeric(args[4])

print(data.frame(indir=indir,outdir=outdir,varj=varj,varj2=varj2))

infiles  <- gtools::mixedsort(list.files(indir,full.names=T))[varj:varj2]
outfiles <- file.path(outdir,basename(infiles))

print("A")
conNames      <- paste0(tools::file_path_sans_ext(basename(infiles)),"_consensus")
print("B")
aln.list      <- lapply(infiles, Biostrings::readDNAMultipleAlignment)
print("C")
aln2.list     <- lapply(aln.list,function(x){Biostrings::DNAMultipleAlignment(gsub("-","N",x))})
print("D")
conStr.list   <- lapply(aln2.list,function(x){Biostrings::consensusString(x, ambiguityMap=IUPAC_CODE_MAP)})
print("E")
conDNA.list   <- lapply(conStr.list,function(x){Biostrings::DNAStringSet(x)})
print("F")
conDNA        <- do.call(c,conDNA.list)
print("G")
names(conDNA) <- conNames
print("F")
out           <- lapply(1:length(outfiles),function(x){Biostrings::writeXStringSet(conDNA[x],outfiles[x])})
