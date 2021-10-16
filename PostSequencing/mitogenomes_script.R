# sbatch --nodes=1 --ntasks-per-node=10 --mem=100Gb --time=6:00:00 --partition=sixhour '/panfs/pfs.local/home/j926w878/work/mitogenomes.sh' <SampleName> <read1 path> [read2 path] [merged reads path]
.libPaths("~/work/R-packages")
library(REEs)

args = commandArgs(trailingOnly=TRUE)
sname      <- args[1]
r1path     <- args[2]

if(args[3]!="NULL"){
	r2path     <- args[3]
} else {
	r2path <- NULL
}
if(args[4]!="NULL"){
	mergedpath <- args[4]
} else {
	mergedpath <- NULL
}

refpath    <- "~/work/mitogenomes/reference.fa"
genespath  <- "~/work/mitogenomes/mtgenes.fa"
outdirpath <- file.path("~/work/mitogenomes/",sname)

bbmappath  <- "~/programs/bbmap/bbmap.sh"
cap3path   <- "~/programs/CAP3/cap3"
spadespath <- "~/programs/SPAdes-3.12.0-Linux/bin/spades.py"
pblatpath  <- "~/work/bi/bin/icebert-pblat-652d3b3/pblat"

result <- get.mitogenome(sampleName=sname,read1=r1path,read2=r2path,reads.merged=mergedpath,referencePATH=refpath,geneFilePATH=genespath,out.dir=outdirpath,bbmapPATH=bbmappath,cap3PATH=cap3path,spadesPATH=spadespath,pblatPATH=pblatpath)
