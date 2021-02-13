# REEs Package

  - [Introduction](#Introduction)
  - [Software Dependencies](#Dependencies)
  - [Installing REEs](#InstallingREEs)
  - [Installing BLAST](#InstallingBLAST)
  - [Installing MAFFT](#InstallingMAFFT)
  - [How to use this package](#HowTo)
  - [Example](#Example)

<a name="Introduction"></a>
## Introduction.
This R package includes functions for selecting **rapidly-evolving exons (REEs)** for targetted sequence capture for phylogenomics. The methods integrated into this package were used to design the [SnakeCap](https://github.com/JeffWeinell/SnakeCap/blob/main/README.md) probe set, a project that was inspired by [FrogCap (Hutter et al.)](https://github.com/chutter/FrogCap-Sequence-Capture) and the [RELEC study of Karin et al. 2019](https://doi.org/10.1093/molbev/msz263).

The reasons for publishing these methods as an R package include (1) having a reproducible and citable pipeline for projects that use the SnakeCap probe set, and (2) to provide a method for researchers to select a set of loci for their phylogenomic studies.

<a name="Dependencies"></a>
### Software Dependencies:
  - [R](https://www.r-project.org/) v3.6+
  - [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - [MAFFT](https://mafft.cbrc.jp/alignment/software/)


**"Quick" install on KU cluster. Delete this before making this repository public.**
```
sbatch --nodes=1 --ntasks-per-node=4 --time=6:00:00 --partition=sixhour /panfs/pfs.local/home/j926w878/programs/REEs_cluster_setup.sh
```


<a name="InstallingREEs"></a>
### Installing REEs
Installation involves the following steps, which are executed in the block of code below: (1–2) define where packages should be installed and where source files should be downloaded, and add these directories to R's library path (3) install and load BiocManager package, which provides access to CRAN and Bioconductor repositories, (4) install dependencies from CRAN and Bioconductor, (5) define URLs to source code of biofiles package (v1.0.0 from CRAN archive) and REEs (SnakeCap project on Open Science Framework), (6) download and install biofiles and REEs from source.

```
# Define where packages should be installed.
# If you want to install or load packages from somewhere other than the default R library path, then set .libPaths() to the directory to use. This must be set each time you start R.

### Probably dont need to use the next few lines. Just set .libPaths().
# packages.dir   <- .libPaths()

# Add packages.dir to the R library path for the duration of the session.
.libPaths("/panfs/pfs.local/scratch/bi/j926w878/scratch_v1/Rv3.6")
#.libPaths("/panfs/pfs.local/scratch/bi/j926w878/scratch_v1/Rv4.0")

# Install BiocManager and several other dependencies
# install.packages(pkgs="BiocManager",lib=packages.dir,repos = "http://cran.us.r-project.org")
#### Test that setting libPaths() works such that lib argument is not needed for install functions, and lib.loc argument not needed for when loading installed packages.
install.packages(pkgs="BiocManager",repos = "http://cran.us.r-project.org")
library(BiocManager)

# Use BiocManager to install REEs and its other dependencies. Packages in their dependency graphs are also installed.
# Set the version argument to "3.10" if using R v3.6; "3.12" for R v4.0; check https://bioconductor.org/about/release-announcements/ to determine which BioConductor version to use for later versions of R.
BiocManager::install(c("BSgenome","DECIPHER","phangorn","dplyr","data.table", "foreach","reutils","curl","knitr","devtools","gschofl/biofiles"),update=FALSE, version="3.12",dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F,Ncpus=4)
# Install REEs. This can be added to the end of the previous line once the REEs repository is public.
BiocManager::install("JeffWeinell/REEs",update=FALSE, version="3.12",dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F,Ncpus=4,auth_token="323d9e4cd00247a39a805dbb66f37db6403cfb8b")
```

<a name="InstallingBLAST"></a>
### Installing BLAST

You can follow the instructions [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/), or use the REEs ```blast.install``` function to install BLAST to the REEs package directory. The ```blast.install``` method is convenient because you wont need to explicitely specify the path to the directory holding BLAST executables when using the REEs ```blast``` function. However, you will need to run ```blast.install``` after each time you upgrade or reinstall REEs.

```
# Reminder: if you closed/reoponed R, remember to set .libPaths() when using a non-default R library location.

# Load REEs
library(REEs)

# Install BLAST
blast.install()
```

<a name="InstallingMAFFT"></a>
### Installing MAFFT

Instructions for installing MAFFT can be found [here](https://mafft.cbrc.jp/alignment/software/).
Also useful: [How to install MAFFT to a non-default location](https://mafft.cbrc.jp/alignment/software/installation_without_root.html).

Alternatively, use the ```install.mafft``` function. The mafft executable will be located "bin" folder of the directory specified by the install.loc argument. If install.loc = "auto" (the default), MAFFT is installed to "/REEs/blast-mafft", and the executable will be in "/REEs/blast-mafft/bin". If install.loc = "PATH", MAFFT is installed to "usr/local", but many users won't have write access to this directory. Keep in mind that if lib.loc="auto" you will need to run ```mafft.install``` after each time that you update or reinstall REEs.
```
# Load REEs
library(REEs)

# Install MAFFT
mafft.install()
```

<a name="HowTo"></a>
## How to use this package to identify rapidly-evolving exons:

Data that you need at the start:
 - A set of genomes (at least four) from the group of interest and preferably from outgroup species. The genomes can be supplied as URLs to genomes on NCBI, or as fasta files.
 - A Generic Feature Format (GFF) table associated with one of the sampled genomes, which can be supplied as a local file or URL path to the file. The species with both a genome and GFF table will be the primary/reference species from which protein-coding sequences will be identified. It is a good idea to use the species with the best annotated and highest quality genome as the reference species. Once you have your genomes and GFF table, work through the pipeline below to select a set of exons to target for sequence capture!

Functions in REEs Pipeline. To see usage details for each function type ```?REEs::<functionName>``` after loading REEs in R.
 1. ```load.gff```
 2. ```filter.gff```
 3. ```get.seqs.from.gff```
 4. ```blast``` (Note: this step is computationally intensive/time consuming)
 5. ```reportBestMatches```
 6. ```get.seqs.from.blastTable```
 7. Repeat steps 4–6 for each genome, with the output of step 3 as the query.
 8. ```makeStatsTable```
 9. ```pick.loci```
 
 The output of step 9 (``pick.loci``) is a ranked list of REEs to target for phylogenomic studies.
 
 The REEs package also provides functions for performing *in-silico* ddRAD, which can be used to asses bias of REEs loci.
 
<a name="Example"></a>
#### Example 

In this example we will identify REEs for the diverse lizard clade [Laterata](https://en.wikipedia.org/wiki/Lacertoidea) (aka Lacertoidea, which includes the Amphisbaenia, Gymnophthalmidae, Lacertidae, and Teiidae). Genomes are available on NCBI for five species of Lacertidae (we will use three of these) and two species of Teiidae. GFF tables are available for three species of Lacertidae: *Lacerta agilis*, *Podarcis muralis*, *Zootoca vivipara*; before we decide which of these three species will be the reference species we will query their GFF tables to get an idea about the quality of the genome and annotation completeness.

Load REEs, define working directory and define the URLs to genomes and GFF tables.
```
### Define R package directory, tempfile directory, and set R library paths variable each time you restart R.
packages.dir <- .libPaths() # Change this if you installed packages somewhere else.
temp.dir     <- tempdir()   # This is probably redundant but doesn't hurt.
.libPaths(c(packages.dir,temp.dir))

### Load REEs package
library(REEs)

### Set path to working directory where files will be saved.
working.dir <- "/panfs/pfs.local/scratch/bi/j926w878/scratch_v1/REEs_example_output"

### URL to GFF of reference species
Podarcis.muralis_GFF.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.gff.gz"
Lacerta.agilis_GFF.url   <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/819/535/GCF_009819535.1_rLacAgi1.pri/GCF_009819535.1_rLacAgi1.pri_genomic.gff.gz"
Zootoca.vivipara_GFF.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/800/845/GCF_011800845.1_UG_Zviv_1/GCF_011800845.1_UG_Zviv_1_genomic.gff.gz"

### URLs to genomes
Podarcis.muralis.genome_url        <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.fna.gz"
Lacerta.agilis.genome_url          <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/819/535/GCF_009819535.1_rLacAgi1.pri/GCF_009819535.1_rLacAgi1.pri_genomic.fna.gz"
Zootoca.vivipara.genome_url        <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/800/845/GCF_011800845.1_UG_Zviv_1/GCF_011800845.1_UG_Zviv_1_genomic.fna.gz"
Aspidoscelis.marmoratus.genome_url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/337/955/GCA_014337955.1_AspMar1.0/GCA_014337955.1_AspMar1.0_genomic.fna.gz"
Salvator.merianae.genome_url       <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/586/115/GCA_003586115.2_HLtupMer6/GCA_003586115.2_HLtupMer6_genomic.fna.gz"
```

<!---
These two genomes are also available but I'm not using them because I've already included a different species of Lacerta.
```
Lacerta.bilineata.genome_url       <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/245/895/GCA_900245895.1_L._bilineata_genome_assembly/GCA_900245895.1_L._bilineata_genome_assembly_genomic.fna.gz"
Lacerta.viridis.genome_url         <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/245/905/GCA_900245905.1_ASM90024590v1/GCA_900245905.1_ASM90024590v1_genomic.fna.gz"
```
--->

Load GFF tables. Count number of scaffolds and number of CDS and genes annotated. Choose which species will be the reference species.
```
### Load GFF
Podarcis.muralis_GFF  <- load.gff(Podarcis.muralis_GFF.url)
Lacerta.agilis_GFF    <- load.gff(Lacerta.agilis_GFF.url)
Zootoca.vivipara_GFF  <- load.gff(Zootoca.vivipara_GFF.url)

### Number of scaffolds
length(unique(Podarcis.muralis_GFF$seqid))  # Result = 2,161
length(unique(Lacerta.agilis_GFF$seqid))    # Result = 29
length(unique(Zootoca.vivipara_GFF$seqid))  # Result = 25,252

# Number of genes annotated
length(grep("gene",Podarcis.muralis_GFF$type))  # Result = 26,745
length(grep("gene",Lacerta.agilis_GFF$type))    # Result = 23,229
length(grep("gene",Zootoca.vivipara_GFF$type))  # Result = 22,898

# Number of CDS annotated
length(grep("CDS",Podarcis.muralis_GFF$type))   # Result = 622,377
length(grep("CDS",Lacerta.agilis_GFF$type))     # Result = 504,650
length(grep("CDS",Zootoca.vivipara_GFF$type))   # Result = 577,940

```
We choose *L. agilis* to be the reference species because sequences are assembled onto fewer scaffolds, likely corresponding to chromosomes. However, *Podarcis muralis* would also be a good choice because it has the most annotate genes and CDS, and sequences are assembled onto relatively few scaffolds. You can work through this pipeline a second time using *P. muralis* as the reference species to compare the REEs selected to those selected using *L. agilis*.

Filter GFF table of the reference species. Extract and save sequences corresponding to regions in the filtered GFF.
<!--- **Note:** I had to modify get.seqs.from.gff because ```contig.names  <- mat.strsplit(headers)[,1]``` will not work if sequence headers have different number of spaces present. I changed this line to ```contig.names  <- mat.strsplit(headers)[,1]```.--->
```
### Filter GFF to include only CDS features ≥ 120bp.
Lacerta.agilis_GFF_CDS_longer120bp <- filter.gff(input.gff=Lacerta.agilis_GFF,feature.type="CDS",min.length=120)

### Save filtered GFF
write.table(x=Lacerta.agilis_GFF_CDS_longer120bp,file=paste0(working.dir,"/Lacerta.agilis_GFF_CDS_longer120bp.txt"),quote=F,sep="\t",row.names=F,col.names=T)

### Extract sequences
Lacerta.agilis_exome   <- get.seqs.from.gff(input.seqs=Lacerta.agilis.genome_url,input.gff=Lacerta.agilis_GFF_CDS_longer120bp)

### Save extracted sequences
Biostrings::writeXStringSet(x=Lacerta.agilis_exome,filepath=paste0(working.dir,"/Lacerta.agilis_exome_longer120bp.fas"))
```

TBLASTX the extracted sequences against each genome. In each case the output is a hit table of matches. The arguments max.matches.per.query sets the maximum number of contigs in the subject species that can contain a match per query sequence; other.args="-max_hsps 10" indicates to return at most 10 matches for each pair of query and subject contigs; eval sets the maximum Expectation Value to consider as a match.
You may need to experiment with the values of these parameters. One strategy is to query a small number of sequences with relaxed settings for -max_hsps, max.matches.per.query, and eval (i.e., use high values for these parameters, e.g. eval=10, max.matches.per.query=50, other.args="-max_hsps 50", to allow many matches per query including low scoring matches). Then examine the output table to determine how stringent (low values) you can set these parameters while still obtaining the likely best match.

For this example, I initially used -max_hsps unlimited, max.matches.per.query=50, and eval=1e-05, but these analyses never finished running and the output files had grown to > 5Gb. Tweaking the parameters to those shown below seemed appropriate, and bitscores were usually high for only the top one to a few matches.
It is a good idea to consult the [BLAST Help Manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/) for a variety of strategies for searching with BLAST. [This table](https://www.ncbi.nlm.nih.gov/books/NBK279684/) of arguments is also useful.

<!---Running with eval=1e-18,max.matches.per.query=7,other.args="-best_hit_overhang 0.1 -best_hit_score_edge 0.1" did not work to keep only the best match per contig pair.--->
```
# Run TBLASTX. I recommend using a bash script and submitting each of these individually to run on a cluster. This will likely take days to run.
Podarcis.muralis.hits        <- REEs::blast(method="tblastx",subject=Podarcis.muralis.genome_url,query=Lacerta.agilis_exome[1:2],table.out=paste0(working.dir,"/Podarcis.muralis.tblastx.exons.hits.txt"),max.targets.per.query=10,max.matches.per.target=10,eval=1e-15)
Lacerta.agilis.hits          <- REEs::blast(method="tblastx",subject=Lacerta.agilis.genome_url,query=Lacerta.agilis_exome,table.out=paste0(working.dir,"/Lacerta.agilis.tblastx.exons.hits.txt"),max.targets.per.query=10,max.matches.per.target=10,eval=1e-15)
Zootoca.vivipara.hits        <- REEs::blast(method="tblastx",subject=Zootoca.vivipara.genome_url,query=Lacerta.agilis_exome,table.out=paste0(working.dir,"/Zootoca.vivipara.tblastx.exons.hits.txt"),max.targets.per.query=10,max.matches.per.target=10,eval=1e-15)
Aspidoscelis.marmoratus.hits <- REEs::blast(method="tblastx",subject=Aspidoscelis.marmoratus.genome_url,query=Lacerta.agilis_exome,table.out=paste0(working.dir,"/Aspidoscelis.marmoratus.tblastx.exons.hits.txt"),max.targets.per.query=10,max.matches.per.target=10,eval=1e-15)
Salvator.merianae.hits       <- REEs::blast(method="tblastx",subject=Salvator.merianae.genome_url,query=Lacerta.agilis_exome,table.out=paste0(working.dir,"/Salvator.merianae.tblastx.exons.hits.txt"),max.targets.per.query=10,max.matches.per.target=10,eval=1e-15)
```

Filter hit tables to remove poor matches (low bitscores) and loci that may have been recently duplicated (i.e., loci having similar bitscores for their best and second best matches, or having at least two matches with bitscores > 60).
```
### Keep the best hit for each match.
best.hits.Podarcis.muralis        <- REEs::reportBestMatches(input.table=Podarcis.muralis.hits,output.table.path=paste0(working.dir,"/Podarcis.muralis.best.hits.txt"))
best.hits.Lacerta.agilis          <- REEs::reportBestMatches(input.table=Lacerta.agilis.hits,output.table.path=paste0(working.dir,"/Lacerta.agilis.best.hits.txt"))
best.hits.Zootoca.vivipara        <- REEs::reportBestMatches(input.table=Zootoca.vivipara.hits,output.table.path=paste0(working.dir,"/Zootoca.vivipara.best.hits.txt"))
best.hits.Aspidoscelis.marmoratus <- REEs::reportBestMatches(input.table=Aspidoscelis.marmoratus.hits,output.table.path=paste0(working.dir,"/Aspidoscelis.marmoratus.best.hits.txt"))
best.hits.Salvator.merianae       <- REEs::reportBestMatches(input.table=Salvator.merianae.hits,output.table.path=paste0(working.dir,"/Salvator.merianae.best.hits.txt"))

```

Get sequences for the best matches from each genome.
```
Podarcis.muralis.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Podarcis.muralis, input.seqs=Podarcis.muralis.genome_url, output.path=paste0(working.dir,"/Podarcis.muralis.tblastx.best.hits_seqs.fas"))
Lacerta.agilis.best.hits.seqs            <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Lacerta.agilis, input.seqs=Lacerta.agilis.genome_url, output.path=paste0(working.dir,"/Lacerta.agilis.tblastx.best.hits_seqs.fas"))
Zootoca.vivipara.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Zootoca.vivipara, input.seqs=Zootoca.vivipara.genome_url, output.path=paste0(working.dir,"/Zootoca.vivipara.tblastx.best.hits_seqs.fas"))
Aspidoscelis.marmoratus.best.hits.seqs   <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Aspidoscelis.marmoratus, input.seqs=Aspidoscelis.marmoratus.genome_url, output.path=paste0(working.dir,"/Aspidoscelis.marmoratus.tblastx.best.hits_seqs.fas"))
Salvator.merianae.best.hits.seqs         <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Salvator.merianae, input.seqs=Salvator.merianae.genome_url, output.path=paste0(working.dir,"/Salvator.merianae.tblastx.best.hits_seqs.fas"))

```

Align homologous sequences and make a table to summarize data in each alignment (one row per aligned locus), and save each alignment to a folder called alignments.
```
### Character vector with paths to "best hits" sequences file of each species.
input.seqs.paths <- paste0(working.dir,c("/Podarcis.muralis.tblastx.best.hits_seqs.fas", "/Lacerta.agilis.tblastx.best.hits_seqs.fas", "/Zootoca.vivipara.tblastx.best.hits_seqs.fas", "/Aspidoscelis.marmoratus.tblastx.best.hits_seqs.fas", "/Salvator.merianae.tblastx.best.hits_seqs.fas"))

### Align sequences and summarize alignments in a table.
stats.table.all  <- makeStatsTable(input.seqs=input.seqs.paths, input.gff=Lacerta.agilis_GFF_CDS_longer120bp, output.path=paste0(working.dir,"/statsTable_REEs.txt"), species= c("Podarcis muralis", "Lacerta agilis", "Zootoca vivipara", "Aspidoscelis marmoratus", "Salvator merianae"),alignments.out=paste0(working.dir,"/alignments"), species.gff=2)
```

Use ```pick.loci``` function to identify REEs to target.
```
### Read the table generated by makeStatsTable function
stats.table.all  <- data.table::fread(paste0(working.dir,"/statsTable_REEs.txt"),header=T)

### Filter the stats table to the optimal set of REEs
stats.table.best <- REEs::pick.loci(statsTable.path=stats.table.all, primary.species="Lacerta agilis", output.path=paste0(working.dir,"/stats_data_FastestExonPerGene_best.tsv"), pident.keep=c(65,100), max.loci.per.gene=1, min.num.species="all", max.capture.coverage=NULL, bait.tiling=0.5, bait.length=120, bait.number=20000, flanking.length="auto", fast.stat="pident", use.min.pident.subgroup=F)
```

**Things to add**
 - Section at beginning of ```pick.loci``` function to check that input arguments have the right class and right combination of input values.
 - Figure out what the plots generate by ```pick.loci``` were intended to show





