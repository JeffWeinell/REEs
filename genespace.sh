module load R
module load anaconda
conda activate py36

### Prepares NCBI GFF3 and protein fasta files into a form that can be read by 

# paths to input GFF and protein fasta files
#GFF_PATH_OG="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/sirtalis_ncbi_dataset/data/GCF_001077635.1/sirtalis_genomic.gff"
#SEQSPATH_OG="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/sirtalis_ncbi_dataset/data/GCF_001077635.1/sirtalis_protein.faa"
GFF_PATH_OG="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/Python-bivittatus_ncbi_datasets/data/GCF_000186305.1/python_genomic.gff"
SEQSPATH_OG="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/Python-bivittatus_ncbi_datasets/data/GCF_000186305.1/python_protein.faa"

# paths to output GFF and protein fasta files
#GFF_PATH_NEW="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/rawGenomes/Tsirtalis/Tsirtalis/annotation/Tsirtalis_gene.gff"
#SEQSPATH_NEW="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/rawGenomes/Tsirtalis/Tsirtalis/annotation/Tsirtalis_pep.fa"
GFF_PATH_NEW="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/rawGenomes/Pbivittatus/Pbivittatus/annotation/Pbivittatus_gene.gff"
SEQSPATH_NEW="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/rawGenomes/Pbivittatus/Pbivittatus/annotation/Pbivittatus_pep.fa"

# create output directory
DATADIR=$(dirname "$SEQSPATH_NEW")
[ ! -d "$DATADIR" ] && mkdir -p "$DATADIR"

#if [ $(awk '{if ($3=="gene"){print}}' "$GFF_PATH_OG" | grep ".*;Name=.*;gene=.*" | wc -l) -gt 1000 ]
#	then
#	KV=$(awk '{if($3=="CDS"){print $9}}' "$GFF_PATH_OG" | awk '{gsub(".*;Name=","Name=")}1' | awk '{gsub(";.*;gene=","\tgene=")}1' | awk '{gsub(";.*","")}1' | grep '.*Name.*gene.*' | sort | uniq | awk '{gsub("Name=","")}1' | awk '{gsub("gene=","")}1' )
#	else
#	KV=$(awk '{if($3=="CDS"){print $9}}' "$GFF_PATH_OG" | awk '{gsub(".*;Name=","Name=")}1' | awk '{gsub(";.*;locus_tag=","\tlocus_tag=")}1' | awk '{gsub(";.*","")}1' | grep '.*Name.*locus_tag.*' | sort | uniq | awk '{gsub("Name=","")}1' | awk '{gsub("locus_tag=","")}1' )
#fi

### KV table for columns for old and new names (use one of the next two lines, depending on the gff)
# KV=$(awk '{if($3=="CDS"){print $9}}' "$GFF_PATH_OG" | awk '{gsub(".*;Name=","Name=")}1' | awk '{gsub(";.*;gene=","\tgene=")}1' | awk '{gsub(";.*","")}1' | grep '.*Name.*gene.*' | sort | uniq | awk '{gsub("Name=","")}1' | awk '{gsub("gene=","")}1' )
# KV=$(awk '{if($3=="CDS"){print $9}}' "$GFF_PATH_OG" | awk '{gsub(".*;Name=","Name=")}1' | awk '{gsub(";.*;locus_tag=","\tlocus_tag=")}1' | awk '{gsub(";.*","")}1' | grep '.*Name.*locus_tag.*' | sort | uniq | awk '{gsub("Name=","")}1' | awk '{gsub("locus_tag=","")}1' )
### This might work for either case above
KV=$(awk '{if($3=="CDS"){print $9}}' "$GFF_PATH_OG" | awk '{gsub(".*;Name=","Name=")}1' | awk '{gsub(";.*;gene=|;.*;locus_tag=","\tgene=")}1' | awk '{gsub(";.*","")}1' | grep '.*Name.*gene.*' | sort | uniq | awk '{gsub("Name=","")}1' | awk '{gsub("gene=","")}1' )

### keep one transcript per gene in protein fasta, and rename sequences as gene name
C1=$(echo "$KV" | awk '{print $1}')
C2=$(echo "$KV" | awk '{print $2}')
KV2=$(Rscript <(echo -e "a = unlist(strsplit(c('"$C1"'),split=' ')) \n b = unlist(strsplit(c('"$C2"'),split=' ')) \n as.data.frame(cbind(a,b)[match(unique(b),b),]) ") | awk 'NR>1{print $2,$3}' | awk '{gsub(" ","\t")}1')
PEP=$(seqkit replace -p '(.+)$' -r '{kv}' -k <(echo "$KV2") <(seqkit replace -p "\s.+" "$SEQSPATH_OG" | seqkit grep -f <(echo "$KV2" | awk '{print $1}')))
echo "$PEP" > "$SEQSPATH_NEW"
### filter the gff to include rows with "gene" in column three, for genes present in the protein fasta; then simplify the gff by removing extra info from column 9
# Simplified GFF
GFF=$(grep -v "^#.*" "$GFF_PATH_OG" | awk '{if($3=="gene"){print}}' | awk '{gsub("ID=.*;gene=|ID=.*;locus_tag=","locus=")}1' | awk '{gsub(";.*","")}1' )
# GFF=$(grep -v "^#.*" "$GFF_PATH_OG" | awk '{if($3=="gene"){print}}' | awk '{gsub("ID=.*;gene=","locus=")}1' | awk '{gsub(";.*","")}1' )
# GFF=$(grep -v "^#.*" "$GFF_PATH_OG" | awk '{if($3=="gene"){print}}' | awk '{gsub("ID=.*;locus_tag=","locus=")}1' | awk '{gsub(";.*","")}1' )
echo "$GFF" > "$GFF_PATH_NEW"

### Remove the original GFF and sequence files
# rm "$SEQSPATH_OG"
# rm "$GFF_PATH_OG"

exit 

################################
module load R
module load anaconda
#conda deactivate py36
conda activate orthofinder
R
.libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
library(GENESPACE)
runwd <-"/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2/"
gpar  <- GENESPACE::init_genespace(
	#genomeIDs=c("Nnaja","Tsirtalis","Telegans","Pbivittatus"),
	#speciesIDs=c("Nnaja","Tsirtalis","Telegans","Pbivittatus"),
	#versionIDs=c("Nnaja","Tsirtalis","Telegans","Pbivittatus"),
	genomeIDs=c("Nnaja","Tsirtalis","Telegans"),
	speciesIDs=c("Nnaja","Tsirtalis","Telegans"),
	versionIDs=c("Nnaja","Tsirtalis","Telegans"),
	ploidy=rep(1,3),
	diamondMode="fast",
	orthofinderMethod="fast",
	wd=runwd,
	nCores=4,
	minPepLen=50,
	gffString="gff",
	pepString="pep",
	path2orthofinder="orthofinder",
	path2mcscanx="/panfs/pfs.local/home/j926w878/programs/MCScanX-master",
	rawGenomeDir=file.path(runwd,"rawGenomes"),
	overwrite=T)

GENESPACE::parse_annotations(
	gsParam = gpar,
	gffEntryType = "gene",
	gffIdColumn = "locus",
	gffStripText = "locus=",
	headerEntryIndex = 1,
	headerSep = " ",
	headerStripText = "locus=",
	overwrite=T)

#gff <- read_gff(gsParam$paths$gff[1])

# number of columns in the gff.gz files; they should all equal 7
ncols_gff <- lengths(lapply(list.files(file.path(runwd,"gff"),full.names=T),function(x){unlist(strsplit(readLines(x,n=1),split="\t"))}))

# gff.gz files with only 6 columns (instead of 7) are missing the 'nbp' column; the statement below adds the nbp column
if (any(ncols_gff==6)){
	gff_paths     <- list.files(file.path(runwd,"gff"),full.names=T)[ncols_gff==6]
	for(i in 1:length(gff_paths)){
		gff_path      <- gff_paths[i]
		gff           <- data.table::fread(gff_path)
		gff[,"nbp"]   <- gff[,"end"]-gff[,"start"]
		gff_path_temp <- gsub(".gz$","",gff_path)
		write.table(gff,file=gff_path_temp,col.names=T,row.names=F,sep="\t",quote=F)
		# remove the old gz gff file
		file.remove(gff_path)
		# gzip the new gff
		system(sprintf("gzip '%s'",gff_path_temp))
	}
}

# run orthofinder to identify orthologous blocks
gpar    <- GENESPACE::run_orthofinder(gsParam = gpar,overwrite=T)
# run synteny to determine synteny among ortholog blocks
gpar2   <- GENESPACE::synteny(gsParam = gpar, overwrite=T)
## 

## plot synteny
ripdat  <- GENESPACE::plot_riparian(gpar2,invertTheseChrs=invert.dt)

#### re-plot synteny with some modifications:
# genome and chromosome names for chomosomes that should be inverted when plotting
invert.dt <- data.table::data.table(genome=c(rep("thamnophis",3),rep("naja",8)),chromosome=c("4","2","5","m2","m3","5","7","m7","m9","m10","m11"))


########
# cp -R "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2_Python/rawGenomes" "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2/rawGenomes"
# rm -R "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/genespace/Naja_Thamnophis2/rawGenomes/Pbivittatus"
########


gffWithOgs.txt.gz
runwd="/Users/jeff/Downloads/"
gpar  <- GENESPACE::init_genespace(
	genomeIDs=c("naja","thamnophis"),
	speciesIDs=c("naja","thamnophis"),
	versionIDs=c("naja","thamnophis"),
	ploidy=rep(1,2),
	diamondMode="fast",
	orthofinderMethod="fast",
	wd=runwd,
	nCores=4,
	minPepLen=50,
	gffString="gff",
	pepString="pep",
	path2orthofinder="orthofinder",
	path2mcscanx=NULL,
	rawGenomeDir=file.path(runwd,"test"))

#########################

gpar=list()
gpar$params$verbose <- T
# gpar$genomes$genomeIDs <- c("naja","thamnophis")
gpar$genomes$genomeIDs <- c("thamnophis","naja")
gpar$genomes$outgroup <- NULL
gpar$paths$results <- "/Users/jeff/Documents/GitHub/SnakeCap"
gpar$params$nCores <- 4

ripdat  <- GENESPACE::plot_riparian(gsParam=gpar, invertTheseChrs=invert.dt, returnSourceData = T)
## The only two files that you need are "gffWithOgs.txt.gz" and "%s_%s_synHits.txt.gz"

ripdat  <- GENESPACE::plot_riparian(gsParam=gpar, invertTheseChrs=invert.dt, returnSourceData = T, plotRegions = F, useOrder = F)



