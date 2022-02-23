########
### sort trees by Naja-naja/Thamnophis-elegans syntenic block
library(ape)
library(IRanges)

#' @param trees multiPhylo object
pRF.dist  <- function(trees){
	grid1 <- expand.grid(1:length(trees),1:length(trees))
	grid2 <- grid1[grid1[,1] < grid1[,2],]
	#ptu   <- apply(grid2,1,function(x){ A=intersect(trees[[x[1]]]$tip.label,trees[[x[2]]]$tip.label); if(length(A)>=4){B=A} else {B=NA}; B}, simplify=FALSE)
	ptu   <- apply(grid2,1,function(x){ intersect(trees[[x[1]]]$tip.label,trees[[x[2]]]$tip.label)}, simplify=FALSE)
	# prf.list <- lapply(1:length(ptu),function(x){ if(any(is.na(ptu[[x]]))) {res=NA} else {t1=trees[[grid2[x,1]]];t2=trees[[grid2[x,2]]];t1b=drop.tip(t1,t1$tip.label[!(t1$tip.label %in% ptu[[x]])]); t2b=drop.tip(t2,t2$tip.label[!(t2$tip.label %in% ptu[[x]])]); res=as.numeric(ape::dist.topo(t1b,t2b))};res})
	prf <- sapply(1:length(ptu),function(x){if(length(ptu[[x]])<4) {res=NA} else {t1=trees[[grid2[x,1]]];t2=trees[[grid2[x,2]]];t1b=drop.tip(t1,t1$tip.label[!(t1$tip.label %in% ptu[[x]])]); t2b=drop.tip(t2,t2$tip.label[!(t2$tip.label %in% ptu[[x]])]); res=as.numeric(ape::dist.topo(unroot(t1b),unroot(t2b)))};res})
	#num.intersect <- unname(lengths(ptu))
	#RF.std <- prf/num.intersect
	if(is.null(names(trees))){
		#names(trees) <- paste0("tree",1:length(trees))
		grid3 <- grid2
	} else {
		grid3 <- names(trees)[as.matrix(grid2)]
		dim(grid3) <- dim(grid2)
	}
	res.df <- data.frame(treeA=grid3[,1],treeB=grid3[,2],RF.int=prf)
	res.df
}


### syntenic blocks were identified with GENESPACE R package
### location of loci in Naja naja was identified with seqkit locate function

# step 1: determine which synteny block each locus belongs to
alncoords.df <- read.table("/Users/jeff/Documents/GitHub/SnakeCap/groups_Naja-naja_genomic_coords_v3.txt",header=F,col.names=c("locus","species","contig","sense","start","end","chr"))
alncoords.df <- alncoords.df[!apply(alncoords.df,1,function(x){any(is.na(x) | x %in% "FAILED" )}),]
mode(alncoords.df[,'start']) <- mode(alncoords.df[,'end']) <- "numeric"
alncoords.df[,'locus'] <- gsub('_v3','',alncoords.df[,'locus'])
alncoords.df[alncoords.df[,'chr']=='na','chr'] <- alncoords.df[alncoords.df[,'chr']=='na','contig']
alncoords.df[,'chr'] <- gsub("MIC_","m",alncoords.df[,'chr'])
synt.df      <- read.table("/Users/jeff/Documents/GitHub/SnakeCap/naja_thamnophis_syntenicBlocks.txt",sep="\t",header=T)
synt.naja.df <- synt.df[synt.df[,"gen1"]=="naja" & synt.df[,"gen2"]=="thamnophis",c("blkID","gen1","chr1","startBp1","endBp1","orient")]
colnames(synt.naja.df) <- c("block","species","chr","start","end","sense")
mode(synt.naja.df[,'start']) <- mode(synt.naja.df[,'end']) <- "numeric"
synt.naja.ldf <- split(synt.naja.df,f= synt.naja.df[,'chr'])
alncoords.ldf <- split(alncoords.df,f= alncoords.df[,'chr'])

synt.naja.ldf2 <- synt.naja.ldf[names(synt.naja.ldf) %in% intersect(names(synt.naja.ldf),names(alncoords.ldf))]
alncoords.ldf2 <- alncoords.ldf[names(alncoords.ldf) %in% intersect(names(synt.naja.ldf),names(alncoords.ldf))]
chromosomes <- sort(names(alncoords.ldf2))

alncoords.ir <- lapply(chromosomes,function(x){IRanges(start=alncoords.ldf2[x][[1]][,'start'], end=alncoords.ldf2[x][[1]][,'end'],names=alncoords.ldf2[x][[1]][,'locus'])})
names(alncoords.ir) <- chromosomes
synt.ir <- lapply(chromosomes,function(x){IRanges(start=synt.naja.ldf2[x][[1]][,'start'], end=synt.naja.ldf2[x][[1]][,'end'],names=synt.naja.ldf2[x][[1]][,'block'])})
names(alncoords.ir) <- chromosomes
######
# finding regions of Naja genome that are not in an orthoblock...
#synt.dir <- lapply(synt.ir,disjoin)
#synt.ir.dir.ol <- lapply(1:length(synt.dir),function(x){findOverlaps(synt.ir[[x]],synt.dir[[x]])})
#test  <- lapply(1:length(synt.dir),function(x){names(synt.ir[[x]])[queryHits(synt.ir.dir.ol[[x]])]})
#test2 <- lapply(test,function(x){grepl("self",x)})
#ol    <- lapply(1:length(chromosomes),function(x){findOverlaps(alncoords.ir[[x]],synt.ir[[x]])})

hits    <- lapply(1:length(synt.ir),function(x){findOverlaps(alncoords.ir[[x]],synt.ir[[x]])})
# checks that each locus is not in more than one orthoblock
all(sapply(1:length(hits),function(x){all(table(queryHits(hits[[x]]))==1)}))

blocknames.new <- lapply(1:length(hits),function(x){paste0("chr",names(alncoords.ir)[x],"block",subjectHits(hits[[x]]))})
blocknames.old <- lapply(1:length(hits),function(x){synt.naja.ldf2[[x]][subjectHits(hits[[x]]),'block']})
locusnames     <- lapply(1:length(hits),function(x){names(alncoords.ir[[x]])[queryHits(hits[[x]])]})
syntenyBlocks_loci <- cbind(unlist(locusnames),unlist(blocknames.new),unlist(blocknames.old))
colnames(syntenyBlocks_loci) <- c("locus","blockname","blockname_genespace")
write.table(syntenyBlocks_loci,"/Users/jeff/Documents/GitHub/SnakeCap/Naja-loci_in_orthoblocks.txt",sep="\t",col.names=T,row.names=F,quote=F)

### summary of results:
# 73 synteny blocks
# 91 loci (of those assigned to a Naja chromosome) were not assigned to any synteny block
# 2262 loci were assigned to a synteny block
###

# step 2: subset gene trees by syntenic block
treespath <- "/Users/jeff/Documents/GitHub/SnakeCap/genetrees/snakecap_3093_genetrees_iqtree_taxanames.trees"
datapath  <- "/Users/jeff/Documents/GitHub/SnakeCap/Naja-loci_in_orthoblocks.txt"
outdir    <- "/Users/jeff/Documents/GitHub/SnakeCap/genetrees/synteny_blocks/"
df        <- read.table(datapath,header=T,sep="\t")
trees.mat <- do.call(rbind,strsplit(readLines(treespath),split=' '))
trees     <- read.tree(text=trees.mat[,2],tree.names=trees.mat[,1])
blocks        <- unique(df$blockname)
trees.list    <- lapply(1:length(blocks),function(x){trees[df[df$blockname %in% blocks[x],'locus']]})
trees.matlist <- lapply(1:length(trees.list),function(x){cbind(names(trees.list[[x]]),write.tree(trees.list[[x]]))})
for(i in 1:length(blocks)){
	outpath <- file.path(outdir,paste0("snakecap_genetrees_iqtree_taxanames_",blocks[i],".trees"))
	write.table(trees.matlist[[i]],file=outpath,sep="\t",col.names=F,row.names=F,quote=F)
}

######
# x=1; test <- cbind(names(trees.list[[x]]),write.tree(trees.list[[x]]))
######

########
# In terminal:

ASTRAL="/Applications/Astral/astral.5.15.5.jar"
OUTDIR="/Users/jeff/Documents/GitHub/SnakeCap/ASTRAL/synteny_blocks"
TREES_INDIR="/Users/jeff/Documents/GitHub/SnakeCap/genetrees/synteny_blocks/"
TREES_INPATHS=$(find $TREES_INDIR -type f | sort -V)

NUMBLOCKS=$(echo "$TREES_INPATHS" | wc -l | awk '{gsub(" ","")}1')
for i in $(seq 1 $NUMBLOCKS ); do

for i in $(seq 11 $NUMBLOCKS); do
  echo $i"/"$NUMBLOCKS
  TREES_INPATH=$(echo "$TREES_INPATHS" | awk -v i="$i" 'NR==i')
  NUMTREES=$(wc -l $TREES_INPATH | awk '{print $1}')
  echo "processing "$NUMTREES" gene trees"
  if [ $NUMTREES -gt 1 ]
  then
    TREE_OG_OUTPATH=$OUTDIR"/astral_"$(basename $TREES_INPATH | awk 'gsub(".*_","")1' | awk 'gsub(".trees",".tree")1' )
    LOG_OG_OUTPATH=$(echo "$TREE_OG_OUTPATH" | awk 'gsub(".tree",".log")1')
    TREES=$(awk '{print $2}' "$TREES_INPATH")
    java -Xmx3000M -jar "$ASTRAL" -i <(echo "$TREES") -o "$TREE_OG_OUTPATH" 2>"$LOG_OG_OUTPATH"
  fi
done

#########
# back in R
library(ape)
astpaths <- list.files("/Users/jeff/Documents/GitHub/SnakeCap/ASTRAL/synteny_blocks/",pattern=".tree",full.names=T)
ast.text <- lapply(astpaths,readLines)
names(ast.text) <- basename(astpaths)
ast.text <- ast.text[lengths(ast.text)>0]
# multiphylo object
ast <- read.tree(text=unlist(ast.text),tree.names=names(ast.text))

### filtering to the set of astral trees containing all tips
tiplabs <- lapply(ast,function(x){x$tip.label})
alllabs <- unique(unlist(tiplabs))
ast.incomplete <- lengths(lapply(tiplabs,function(x){setdiff(alllabs,x)}))>0
ast2 <- ast[!ast.incomplete]
ast3 <- lapply(ast2,function(x){root(x,node=MRCA(x,c("Python-molurus","Boa-constrictor")))})
ast4 <- read.tree(text=sapply(ast3,write.tree),tree.names=names(ast3))

### calculating RF distances between astral trees
rf      <- dist.topo(ast4)
rf.list <- as.list(as.data.frame(as.matrix(rf)))
rf.df   <- do.call(rbind,lapply(1:length(rf.list),function(x){data.frame(tree1=names(rf.list),tree2=names(rf.list)[x],RF=unname(rf.list[[x]]))}))
rf.df2  <- rf.df[rf.df[,'tree1']<rf.df[,'tree2'],]
write.table(rf.df2, "/Users/jeff/Documents/GitHub/SnakeCap/RF_distances_ASTRAL_synenty-blocks.txt",col.names=T,sep="\t",quote=F,row.names=F)

###
# results: none of the within-synteny block astral trees have the same topology
###

rf.df2[,"chr.tree1"] <- gsub("astral_chr|block.*","",rf.df2[,"tree1"])
rf.df2[,"chr.tree2"] <- gsub("astral_chr|block.*","",rf.df2[,"tree2"])
# rf.df2[,"chr.tree1"] <- gsub("astral_chr|block.*","",rf.df2[,"tree1"])

###################################
# ASTRAL tree for each chromosome #
# subset gene trees by chromosome
treespath <- "/Users/jeff/Documents/GitHub/SnakeCap/genetrees/snakecap_3093_genetrees_iqtree_taxanames.trees"
outdir    <- "/Users/jeff/Documents/GitHub/SnakeCap/genetrees/chromosomes_naja/"
datapath  <- "/Users/jeff/Documents/GitHub/SnakeCap/groups_Naja-naja_genomic_coords_v3.txt"
df        <- read.table(datapath,header=F,col.names=c("locus","species","contig","sense","start","end","chr"))
df <- df[!apply(df,1,function(x){any(is.na(x) | x %in% "FAILED" )}),]
mode(df[,'start']) <- mode(df[,'end']) <- "numeric"
df[,'locus'] <- gsub('_v3','',df[,'locus'])
df <- df[df[,'chr']!='na',]
df[,'chr'] <- gsub("MIC_","m",df[,'chr'])

trees.mat <- do.call(rbind,strsplit(readLines(treespath),split=' '))
trees     <- read.tree(text=trees.mat[,2],tree.names=trees.mat[,1])
chrom <- unique(df$chr)
trees.list    <- lapply(1:length(chrom),function(x){trees[df[df$chr %in% chrom[x],'locus']]})
trees.matlist <- lapply(1:length(trees.list),function(x){cbind(names(trees.list[[x]]),write.tree(trees.list[[x]]))})
for(i in 1:length(chrom)){
	outpath <- file.path(outdir,paste0("snakecap_genetrees_iqtree_taxanames_chromosome-",chrom[i],".trees"))
	write.table(trees.matlist[[i]],file=outpath,sep="\t",col.names=F,row.names=F,quote=F)
}

########
# In terminal:

ASTRAL="/Applications/Astral/astral.5.15.5.jar"
OUTDIR="/Users/jeff/Documents/GitHub/SnakeCap/ASTRAL/chromosomes_naja"
TREES_INDIR="/Users/jeff/Documents/GitHub/SnakeCap/genetrees/chromosomes_naja/"
TREES_INPATHS=$(find $TREES_INDIR -type f | sort -V)

NUMCHROM=$(echo "$TREES_INPATHS" | wc -l | awk '{gsub(" ","")}1')
for i in $(seq 2 $NUMCHROM ); do
  echo $i"/"$NUMCHROM
  TREES_INPATH=$(echo "$TREES_INPATHS" | awk -v i="$i" 'NR==i')
  NUMTREES=$(wc -l $TREES_INPATH | awk '{print $1}')
  echo "processing "$NUMTREES" gene trees"
  if [ $NUMTREES -gt 1 ]
  then
    TREE_OG_OUTPATH=$OUTDIR"/astral_"$(basename $TREES_INPATH | awk 'gsub(".*_","")1' | awk 'gsub(".trees",".tree")1' )
    LOG_OG_OUTPATH=$(echo "$TREE_OG_OUTPATH" | awk 'gsub(".tree",".log")1')
    TREES=$(awk '{print $2}' "$TREES_INPATH")
    java -Xmx3000M -jar "$ASTRAL" -i <(echo "$TREES") -o "$TREE_OG_OUTPATH" 2>"$LOG_OG_OUTPATH"
  fi
done
#######
# Back in R:
library(ape)
astpaths <- list.files("/Users/jeff/Documents/GitHub/SnakeCap/ASTRAL/chromosomes_naja/",pattern=".tree",full.names=T)
ast.text <- lapply(astpaths,readLines)
names(ast.text) <- basename(astpaths)
ast.text <- ast.text[lengths(ast.text)>0]
# multiphylo object
ast <- read.tree(text=unlist(ast.text),tree.names=names(ast.text))

### filtering to the set of astral trees containing all tips
tiplabs <- lapply(ast,function(x){x$tip.label})
alllabs <- unique(unlist(tiplabs))
ast.incomplete <- lengths(lapply(tiplabs,function(x){setdiff(alllabs,x)}))>0
ast2 <- ast[!ast.incomplete]
ast3 <- lapply(ast2,function(x){root(x,node=MRCA(x,c("Python-molurus","Boa-constrictor")))})
ast4 <- read.tree(text=sapply(ast3,write.tree),tree.names=names(ast3))

### calculating RF distances between astral trees
rf      <- dist.topo(ast4)
rf.list <- as.list(as.data.frame(as.matrix(rf)))
rf.df   <- do.call(rbind,lapply(1:length(rf.list),function(x){data.frame(tree1=names(rf.list),tree2=names(rf.list)[x],RF=unname(rf.list[[x]]))}))
rf.df2  <- rf.df[rf.df[,'tree1']<rf.df[,'tree2'],]
write.table(rf.df2, "/Users/jeff/Documents/GitHub/SnakeCap/RF_distances_ASTRAL_chromosomes_naja.txt",col.names=T,sep="\t",quote=F,row.names=F)

### heat matrix for RF distances among chromosomes?

# testing the effect of dropping Buhoma... still no chromosome trees with RF dist = 0
ast5    <- drop.tip.multiPhylo(ast4,"Buhoma-depressiceps_CFS1528")
rf      <- dist.topo(ast5)
rf.list <- as.list(as.data.frame(as.matrix(rf)))
rf.df   <- do.call(rbind,lapply(1:length(rf.list),function(x){data.frame(tree1=names(rf.list),tree2=names(rf.list)[x],RF=unname(rf.list[[x]]))}))
rf.df2  <- rf.df[rf.df[,'tree1']<rf.df[,'tree2'],]

for(i in 1:length(ast5)){
	ast5[[i]]$edge.length <- rep(1,length(ast5[[i]]$edge.length))
}

# ggtree(ast5["astral_chromosome-1.tree"]) + geom_tiplab(size=2) + xlim(0,30)

##########################
## Binning by position along chromosome:
chrpath <- "/Users/jeff/Documents/GitHub/SnakeCap/naja_chromosomes.txt" 
chr.df <- read.table(chrpath,sep="\t",header=F,col.names=c("chr","start","end"))

ranges.outer <- data.frame(start=1,ceiling(chr.df[,'end']/4)

start1 <- chr.df[,'start']
end1   <- ceiling(chr.df[,'end']/4)
start2 <- end1+1
end2   <- end1*3
start3 <- end2+1
end3   <- chr.df[,'end']

chr.df[,'chr'] <- gsub("MIC_","m",chr.df[,'chr'])
ranges1 <- IRanges(start1,end1,names=chr.df[,'chr'])
ranges2 <- IRanges(start2,end2,names=chr.df[,'chr'])
ranges3 <- IRanges(start3,end3,names=chr.df[,'chr'])

hits1 <- lapply(1:length(chromosomes),function(x){i1=which(names(alncoords.ir) %in% chromosomes[x]);i2=which(names(ranges1) %in% chromosomes[x]); findOverlaps(alncoords.ir[[i1]],ranges1[i2])})
hits2 <- lapply(1:length(chromosomes),function(x){i1=which(names(alncoords.ir) %in% chromosomes[x]);i2=which(names(ranges2) %in% chromosomes[x]); findOverlaps(alncoords.ir[[i1]],ranges2[i2])})
hits3 <- lapply(1:length(chromosomes),function(x){i1=which(names(alncoords.ir) %in% chromosomes[x]);i2=which(names(ranges3) %in% chromosomes[x]); findOverlaps(alncoords.ir[[i1]],ranges3[i2])})

hits1_loci <- lapply(1:length(hits1),function(x){names(alncoords.ir[[x]])[queryHits(hits1[[x]])]})
hits2_loci <- lapply(1:length(hits2),function(x){names(alncoords.ir[[x]])[queryHits(hits2[[x]])]})
hits3_loci <- lapply(1:length(hits3),function(x){names(alncoords.ir[[x]])[queryHits(hits3[[x]])]})

outer.loci <- lapply(1:length(hits1_loci),function(x){c(hits1_loci[[x]],hits3_loci[[x]])})
inner.loci <- hits2_loci

outer.df <- lapply(1:length(alncoords.ir),function(x){data.frame(locus=outer.loci[[x]],chr=names(alncoords.ir)[x],location="outer")})
inner.df <- lapply(1:length(alncoords.ir),function(x){data.frame(locus=inner.loci[[x]],chr=names(alncoords.ir)[x],location="inner")})
outer.df <- do.call(rbind, outer.df)
inner.df <- do.call(rbind, outer.df)
locs.df <- rbind(outer.df,inner.df)
locs.df[,'location'] <- paste0("chr",locs.df[,'chr'],"_",locs.df[,'location'])
write.table(locs.df,"/Users/jeff/Documents/GitHub/SnakeCap/loci_chromosome-locations_bigBins.txt",sep="\t",col.names=T,row.names=F,quote=F)

#####################################################################
#### sort trees into outer and inner groups for each chromosome #####
#####################################################################
treespath <- "/Users/jeff/Documents/GitHub/SnakeCap/genetrees/snakecap_3093_genetrees_iqtree_taxanames.trees"
outdir    <- "/Users/jeff/Documents/GitHub/SnakeCap/genetrees/chromosomes_outer_vs_inner/"
datapath  <- "/Users/jeff/Documents/GitHub/SnakeCap/loci_chromosome-locations_bigBins.txt"
df        <- read.table(datapath,header=T,sep="\t")

trees.mat <- do.call(rbind,strsplit(readLines(treespath),split=' '))
trees     <- read.tree(text=trees.mat[,2],tree.names=trees.mat[,1])
loc       <- gtools::mixedsort(unique(df$location))
trees.list    <- lapply(1:length(loc),function(x){trees[df[df$location %in% loc[x],'locus']]})
trees.matlist <- lapply(1:length(trees.list),function(x){cbind(names(trees.list[[x]]),write.tree(trees.list[[x]]))})
for(i in 1:length(loc)){
	outpath <- file.path(outdir,paste0("snakecap_genetrees_iqtree_taxanames_",loc[i],".trees"))
	write.table(trees.matlist[[i]],file=outpath,sep="\t",col.names=F,row.names=F,quote=F)
}

pRFdist.list <- lapply(trees.list,pRF.dist)
names(pRFdist.list) <- loc
RF.loc.means <- sapply(pRFdist.list,function(x){mean(x[,'RF.int'],na.rm=T)})

########
# In terminal:
ASTRAL="/Applications/Astral/astral.5.15.5.jar"
OUTDIR="/Users/jeff/Documents/GitHub/SnakeCap/ASTRAL/chromosomes_outer_vs_inner"
TREES_INDIR="/Users/jeff/Documents/GitHub/SnakeCap/genetrees/chromosomes_outer_vs_inner/"
TREES_INPATHS=$(find $TREES_INDIR -type f | sort -V)

NUM=$(echo "$TREES_INPATHS" | wc -l | awk '{gsub(" ","")}1')
for i in $(seq 1 $NUM ); do
  echo $i"/"$NUM
  TREES_INPATH=$(echo "$TREES_INPATHS" | awk -v i="$i" 'NR==i')
  NUMTREES=$(wc -l $TREES_INPATH | awk '{print $1}')
  echo "processing "$NUMTREES" gene trees"
  if [ $NUMTREES -gt 1 ]
  then
    TREE_OG_OUTPATH=$OUTDIR"/astral_"$(basename $TREES_INPATH | awk 'gsub(".*taxanames_","")1' | awk 'gsub(".trees",".tree")1' )
    LOG_OG_OUTPATH=$(echo "$TREE_OG_OUTPATH" | awk 'gsub(".tree",".log")1')
    TREES=$(awk '{print $2}' "$TREES_INPATH")
    java -Xmx3000M -jar "$ASTRAL" -i <(echo "$TREES") -o "$TREE_OG_OUTPATH" 2>"$LOG_OG_OUTPATH"
  fi
done
#######
# In R:
astpaths <- list.files("/Users/jeff/Documents/GitHub/SnakeCap/ASTRAL/chromosomes_outer_vs_inner/",pattern=".tree",full.names=T)
ast.text <- lapply(astpaths,readLines)
names(ast.text) <- basename(astpaths)
ast.text <- ast.text[lengths(ast.text)>0]
# multiphylo object
ast <- read.tree(text=unlist(ast.text),tree.names=names(ast.text))

### filtering to the set of astral trees containing all tips
tiplabs <- lapply(ast,function(x){x$tip.label})
alllabs <- unique(unlist(tiplabs))
ast.incomplete <- lengths(lapply(tiplabs,function(x){setdiff(alllabs,x)}))>0
ast2 <- ast[!ast.incomplete]
ast3 <- lapply(ast2,function(x){root(x,node=MRCA(x,c("Python-molurus","Boa-constrictor")))})
ast4 <- read.tree(text=sapply(ast3,write.tree),tree.names=names(ast3))

### calculating RF distances between astral trees
rf      <- dist.topo(ast4)
rf.list <- as.list(as.data.frame(as.matrix(rf)))
rf.df   <- do.call(rbind,lapply(1:length(rf.list),function(x){data.frame(tree1=names(rf.list),tree2=names(rf.list)[x],RF=unname(rf.list[[x]]))}))
rf.df2  <- rf.df[rf.df[,'tree1']<rf.df[,'tree2'],]
# write.table(rf.df2, "/Users/jeff/Documents/GitHub/SnakeCap/RF_distances_ASTRAL_chromosomes_naja.txt",col.names=T,sep="\t",quote=F,row.names=F)

####################################
datapath  <- "/Users/jeff/Documents/GitHub/SnakeCap/groups_Naja-naja_genomic_coords_v3.txt"
df        <- read.table(datapath,header=F,col.names=c("locus","species","contig","sense","start","end","chr"))
df <- df[!apply(df,1,function(x){any(is.na(x) | x %in% "FAILED" )}),]
mode(df[,'start']) <- mode(df[,'end']) <- "numeric"
df[,'locus'] <- gsub('_v3','',df[,'locus'])
df <- df[df[,'chr']!='na',]
df[,'chr'] <- gsub("MIC_","m",df[,'chr'])
chromosomes     <- gtools::mixedsort(unique(unlist(df[,'chr'])))
chromosomes.key <- 1:length(chromosomes)
df[,'chr.id2'] <- chromosomes.key[match(df[,'chr'], chromosomes)]
df2            <- df[order(df[,'chr.id2'],df[,'start']),]
chr.pos <- unname(unlist(lapply(table(df2[,'chr.id2']),function(x){1:x})))
df2[,'locus.id2'] <- paste0("C",df2[,'chr.id2'],"L",chr.pos)
### only includes loci that could be mapped to Naja chromosomes... so 2311 loci
write.table(df2,"/Users/jeff/Documents/GitHub/SnakeCap/loci_coords_ordered_naja-chromosmes.txt",col.names=T,row.names=F,sep="\t",quote=F)
####################################





