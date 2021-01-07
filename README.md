# SnakeCap Sequence Capture Probe Set

# Contents

[Description of SnakeCap Probe Set](#Description)

[Methods](#Methods)
  - [Choosing loci to capture](#Methods.SelectingTargetLoci)
    - [Rapidly-evolving Exons (REEs)](#Methods.SelectingREEs)
    - [Ultraconserved-elements (UCEs)](#Methods.SelectingUCEs)
    - [ddRAD-like loci](#Methods.SelectingddRAD)
    - [major histocompatibility loci (MHCs)](#Methods.SelectingMHC)
    - [scalation loci](#Methods.SelectingScalation)
    - [vision loci](#Methods.SelectingVision)
  - [Probe Synthesis](#ProbeSynthesis)
  - [Taxa Sampled](#Sampling)
  - [Sequence Capture Library Prep](#LibraryPrep)
  - [DNA Sequencing](#DNASequencing)
  - [Post-sequencing](#PostSequencing)
    - [Demultiplexing](#Demultiplexing)
    - [Processing sequence reads](#ProcessingReads)
    - [DNA Alignment](#DNA.Alignment)
    - [Phylogenetic Analyses](#PhylogeneticAnalyses)

[Results](#Results)

[References](#References)

<a name="Description"></a>
# Description of SnakeCap Probe (Bait) Set

The probe set includes 20,020 probes for 3,129 single-copy loci (1,517,011 nt) shared across snakes. The target loci are categorized into four types: (1) rapidly evolving exons (REEs; n = 1,653), (2) ultra-conserved elements (UCEs; n = 907), (3) ddRAD-like loci (n = 328), (4) and functionally interesting genes, which includes 27 major histocompatibility complex (MHC) genes, 119 vision genes, and 95 scalation genes.

REEs include one or more entire exons and one or both exon-flanking regions, and range in length from 121 to 7,501 nt. I used a modified version of the FrogCap pipeline (Hutter et al., 2019) to select the optimal set of REEs from an alignment of snake exomes.

SnakeCap UCEs are a subset of the *Micrurus fulvius* UCEs from Streicher and Wiens (2017).

ddRAD-like loci are shared, single-copy loci identified from in-silico ddRAD using recognition sites for SbfI and EcoRI restriction enzymes.

Functional loci included entire or partial gene regions that have previously been predicted or known to function in either (1) vertebrate immune systems, (2) vision, (3) or scalation.

#### Table 1. For each type of locus: genomic region targeted, number of loci (nloci), total number of nucleotides targeted (nnt), and nucleotide lengths (nt/locus) of the shortest and longest loci.

Locus type | Region targeted | nloci | nnt | nt/locus (min–max)
---- | ---- | ---- | ---- | ----
Rapidly evolving exons (REEs) | Usually the entire exon + 0–60nt each of 5' & 3' flanking regions. | 1,653 | 996,369 | 121–7,501
Ultra-conserved elements (UCEs) | Entire UCE region previously identified | 907 | 143,075 | 120–161
ddRAD-like | in silico ddRAD selected loci | 328 | 271,505 | 120–996
MHC | Exons + 0–60nt each of 5' & 3' flanking regions of major histocompatibility complex (MHC) genes. | 27 | 5,354 | 121–364
Vision | 160nt, including ≤ 70nt upstream of start codon if targeted exon is first exon. Usually, first 160nt of exon targeted. | 119 | 18,857 | 120–170
Scalation | 1100nt, including 1000nt of promoter region + first 100nt of first exon. | 95 | 81,851 | 125–1,101
All loci |  | 3,129 | 1,517,011 | 120–7,501 (mean = 531.62)  


<a name="Methods"></a>
# Methods

<a name="Methods.SelectingREEs"></a>
## Selecting the set of target REEs

<!-- <a name="Methods.SelectingREEs.overview"></a> -->
#### Overview: 

I used the following R functions (shown as packageName::functionName):

REEs::load.gff --> REEs::filter.gff --> REEs::get.seqs.from.gff --> Biostrings::writeXStringSet --> REEs::blast --> REEs::reportBestMatches --> REEs::get.seqs.from.blastTable --> REEs::makeStatsTable --> REEs::pick.loci --> (then functions in step 8  to get REEs + small region of exon-flanking DNA).

Further description and implementation of these functions is shown below in Details section.

<!-- <a name="Methods.SelectingREEs.detailed"></a> -->
#### Details:

<!--
I downloaded the [*Thamnophis sirtalis* genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz) and its associated annotation table: [GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz) (n = 559,130 features annotated). Then, I renamed the contigs in the genome file to have the following format: **Thamnophis_sirtalis_GCF_001077635.1_read1**, **Thamnophis_sirtalis_GCF_001077635.1_read2**, etc., and saved this renamed genome in sequential fasta format: [ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3.zip](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/exomes/ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3.zip). The two-column, tab-delimited table [Scaffold-Name-Key.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/exomes/Scaffold-Name-Key.txt?token=AJJOG2UQ6MDA7UY2U4R6BFS7ZDZYS) includes the new contig name in the first column and the original contig name in the second column:
-->

1. I used the functions load.gff and filter.gff (REEs R package) to filter the *Thamnophis sirtalis* genome feature table (GFF3 format) to include only CDS feature tracks ≥ 120bp in length. The input (unfiltered) GFF3 file can be found here: [GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz). Learn about GFF3 file format [here](https://uswest.ensembl.org/info/website/upload/gff3.html).

```
### Load REEs package
library(REEs)

### Define URL to the Thamnophis sirtalis genome feature table.
Thamnophis.sirtalis_GFF.url       <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"

### Load the input feature table into R
Thamnophis.sirtalis_GFF     <- load.gff(input=Thamnophis.sirtalis_GFF.url,local=F)

### Filter Thamnophis.sirtalis_GFF to only include CDS features with length at least 120bp (the size of the baits).
Thamnophis.sirtalis_GFF_CDS_longer120bp <- filter.gff(input.gff=Thamnophis.sirtalis_GFF,feature.type="CDS",min.length=120)

### Save the filtered feature table
write.table(x=output.gff,file="./Thamnophis.sirtalis_GFF_CDS_longer120bp.txt",quote=F,sep="\t",row.names=F,col.names=T)
```
The output feature table can be downloaded here: [Thamnophis.sirtalis_GFF_CDS_longer120bp.txt](https://osf.io/tm7qw/download).  Note 1: This filtered GFF table is also included as a data table object in the REEs R package (object name: Thamnophis.sirtalis_GFF_CDS_longer120bp). Note 2: the filtered GFF table is not true GFF format, because the header/comment lines are not included; nevertheless, the format used for the columns follows GFF3 format.

2. I used the function get.seqs.from.gff to extract the sequences corresponding the CDS features in the table generated in step 1. 

```
### URL to the Thamnophis sirtalis genome (fasta formatted sequences).
Thamnophis.sirtalis_genome.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz"

### Extract the sequences for the loci in the filtered GFF (Thamnophis.sirtalis_GFF_CDS_longer120bp from step 1).
Thamnophis.sirtalis_exome   <- get.seqs.from.gff(input.seqs=Thamnophis.sirtalis_genome.path,input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp) 

### Save extracted sequences
writeXStringSet(x=Thamnophis.sirtalis_exome,filepath="./Thamnophis_sirtalis_exome_longer120bp.fas")

```
Output DNA sequences were saved in fasta format and can be downloaded here: [Thamnophis_sirtalis_exome_longer120bp.fas](https://osf.io/v7ecb/download).

3. I used TBLASTX as implemented in the REEs::blast function to search for each *T. sirtalis* CDS sequence (from step 2) in each squamate genome listed in Table 2. I saved up to 50 matches per query sequence.

**Table 2**. Genomes used to select REEs included all squamate genomes available from NCBI in 2017.
Species  | Family | NCBI Genome Assembly Accession
----|----|---- 
*Anolis carolinensis* | Dactyloidae | GCF_000090745.1
*Gekko japonicus* | Gekkonidae | GCF_001447785.1
*Pogona vitticeps* | Agamidae | GCF_900067755.1
*Crotalus horridus* | Viperidae (Crotalinae) | GCA_001625485.1
*Crotalus mitchellii* (= *C. pyrrhus*) | Viperidae (Crotalinae) | GCA_000737285.1
*Ophiophagus hannah* | Elapidae | GCA_000516915.1
*Pantherophis guttatus* | Colubridae (Colubrinae) | GCA_001185365.1
*Protobothrops mucrosquamatus* | Viperidae (Crotalinae) | GCF_001527695.1
*Python bivitattus* | Pythonidae | GCF_000186305.1
*Thamnophis sirtalis* | Colubridae (Natricinae) | GCF_001077635.1
*Vipera berus berus* | Viperidae (Viperinae) | GCA_000800605.1

```
### URLs to the genomes used are held in a matrix that can be accessed with the REEs::datasets function.
Anolis.carolinensis.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Anolis carolinensis"),2]
Gekko.japonicus.genome_url              <- REEs::datasets(1)[which(datasets(1)[,1]=="Gekko japonicus"),2]
Pogona.vitticeps.genome_url             <- REEs::datasets(1)[which(datasets(1)[,1]=="Pogona vitticeps"),2]
Crotalus.horridus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus horridus"),2]
Crotalus.mitchellii.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus mitchellii"),2]
Ophiophagus.hannah.genome_url           <- REEs::datasets(1)[which(datasets(1)[,1]=="Ophiophagus hannah"),2]
Pantherophis.guttatus.genome_url        <- REEs::datasets(1)[which(datasets(1)[,1]=="Pantherophis guttatus"),2]
Protobothrops.mucrosquamatus.genome_url <- REEs::datasets(1)[which(datasets(1)[,1]=="Protobothrops mucrosquamatus"),2]
Python.bivittatus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Python bivittatus"),2]
Viperus.berus.genome_url                <- REEs::datasets(1)[which(datasets(1)[,1]=="Viperus berus"),2]
Thamnophis.sirtalis.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Thamnophis sirtalis"),2]

# Runs TBLASTX
Anolis.carolinensis.50hits          <- REEs::blast(method="tblastx",subject=Anolis.carolinensis.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Anolis.carolinensis.tblastx.exons.50hits.txt")
Gekko.japonicus.50hits              <- REEs::blast(method="tblastx",subject=Gekko.japonicus.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Gekko.japonicus.tblastx.exons.50hits.txt")
Pogona.vitticeps.50hits             <- REEs::blast(method="tblastx",subject=Pogona.vitticeps.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Pogona.vitticeps.tblastx.exons.50hits.txt")
Crotalus.horridus.50hits            <- REEs::blast(method="tblastx",subject=Crotalus.horridus.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Crotalus.horridus.tblastx.exons.50hits.txt")
Crotalus.mitchellii.50hits          <- REEs::blast(method="tblastx",subject=Crotalus.mitchellii.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Crotalus.mitchellii.tblastx.exons.50hits.txt")
Ophiophagus.hannah.50hits           <- REEs::blast(method="tblastx",subject=Ophiophagus.hannah.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Ophiophagus.hannah.tblastx.exons.50hits.txt")
Pantherophis.guttatus.50hits        <- REEs::blast(method="tblastx",subject=Pantherophis.guttatus.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Pantherophis.guttatus.tblastx.exons.50hits.txt")
Protobothrops.mucrosquamatus.50hits <- REEs::blast(method="tblastx",subject=Protobothrops.mucrosquamatus.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Protobothrops.mucrosquamatus.tblastx.exons.50hits.txt")
Python.bivittatus.50hits            <- REEs::blast(method="tblastx",subject=Python.bivittatus.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Python.bivittatus.tblastx.exons.50hits.txt")
Viperus.berus.50hits                <- REEs::blast(method="tblastx",subject=Viperus.berus.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Viperus.berus.tblastx.exons.50hits.txt")
Thamnophis.sirtalis.50hits          <- REEs::blast(method="tblastx",subject=Thamnophis.sirtalis.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Thamnophis.sirtalis.tblastx.exons.50hits.txt")
```
The TBLASTX output tables can be downloaded here (these are xz compressed): [Anolis.carolinensis.tblastx.exons.50hits.txt.xz](https://osf.io/3a4be/download), [Gekko.japonicus.tblastx.exons.50hits.txt.xz](https://osf.io/u2ec5/download), [Pogona.vitticeps.tblastx.exons.50hits.txt.xz](https://osf.io/9ehyf/download), [Crotalus.horridus.tblastx.exons.50hits.txt.xz](https://osf.io/j4gvk/download), [Crotalus.mitchellii.tblastx.exons.50hits.txt.xz](https://osf.io/5vxt7/download), [Ophiophagus.hannah.tblastx.exons.50hits.txt.xz](https://osf.io/jr5zd/download), [Pantherophis.guttatus.tblastx.exons.50hits.txt.xz](https://osf.io/tu3ve/download), [Protobothrops.mucrosquamatus.tblastx.exons.50hits.txt.xz](https://osf.io/7xmr4/download), [Python.bivittatus.tblastx.exons.50hits.txt.xz](https://osf.io/8s5na/download), [Viperus.berus.tblastx.exons.50hits.txt.xz](https://osf.io/48x7e/download), [Thamnophis.sirtalis.tblastx.exons.50hits.txt.xz](https://osf.io/ctz73/download).

4. I used the REEs::reportBestMatches function to filter the hit tables generated by TBLASTX to include only the best match (max bitscore) per query. Matches with bitscore < 50 were dropped. If bitscores of the best and second-best matches differed by < **5**, then the query locus was filtered.

```
best.hits.Anolis.carolinensis           <- REEs::reportBestMatches(input.table=Anolis.carolinensis.50hits,output.table.path="./Anolis.carolinensis.best.hits.txt")
best.hits.Gekko.japonicus               <- REEs::reportBestMatches(input.table=Gekko.japonicus.50hits,output.table.path="./Gekko.japonicus.best.hits.txt")
best.hits.Pogona.vitticeps              <- REEs::reportBestMatches(input.table=Pogona.vitticeps.50hits,output.table.path="./Pogona.vitticeps.best.hits.txt")
best.hits.Crotalus.horridus             <- REEs::reportBestMatches(input.table=Crotalus.horridus.50hits,output.table.path="./Crotalus.horridus.best.hits.txt")
best.hits.Crotalus.mitchellii           <- REEs::reportBestMatches(input.table=Crotalus.mitchellii.50hits,output.table.path="./Crotalus.mitchellii.best.hits.txt")
best.hits.Ophiophagus.hannah            <- REEs::reportBestMatches(input.table=Ophiophagus.hannah.50hits,output.table.path="./Ophiophagus.hannah.best.hits.txt")
best.hits.Pantherophis.guttatus         <- REEs::reportBestMatches(input.table=Pantherophis.guttatus.50hits,output.table.path="./Pantherophis.guttatus.best.hits.txt")
best.hits.Protobothrops.mucrosquamatus  <- REEs::reportBestMatches(input.table=Protobothrops.mucrosquamatus.50hits,output.table.path="./Protobothrops.mucrosquamatus.best.hits.txt")
best.hits.Python.bivittatus             <- REEs::reportBestMatches(input.table=Python.bivittatus.50hits,output.table.path="./Python.bivittatus.best.hits.txt")
best.hits.Viperus.berus                 <- REEs::reportBestMatches(input.table=Viperus.berus.50hits,output.table.path="./Viperus.berus.best.hits.txt")
```

Output tables containing the set of best matches: [Anolis.carolinensis.exons.best.hits.txt](https://osf.io/a5b8n/download), [Gekko.japonicus.exons.best.hits.txt](https://osf.io/a9p23/download), [Pogona.vitticeps.exons.best.hits.txt](https://osf.io/xrycw/download), [Crotalus.horridus.exons.best.hits.txt](https://osf.io/zrndt/download), [Crotalus.mitchellii.exons.best.hits.txt](https://osf.io/bhnwy/download), [Ophiophagus.hannah.exons.best.hits.txt](https://osf.io/y49cv/download), [Pantherophis.guttatus.exons.best.hits.txt](https://osf.io/2vw9f/download), [Protobothrops.mucrosquamatus.exons.best.hits.txt](https://osf.io/aqps2/download), [Python.bivittatus.exons.best.hits.txt](https://osf.io/wxthd/download), [Viperus.berus.exons.best.hits.txt](https://osf.io/en4yg/download),[Thamnophis.sirtalis.exons.best.hits.txt](https://osf.io/mekqz/download).

5. I used get.seqs.from.blastTable function (REEs package) to extract sequences for the best matches identified in step 4.

```
#### Extracts the sequence of the best match of each exon query from each genome.
Anolis.carolinensis.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Anolis.carolinensis,input.seqs=Anolis.carolinensis.genome_url,output.path="./Anolis.carolinensis.tblastx.best.hits_seqs.fas")
Gekko.japonicus.best.hits.seqs              <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Gekko.japonicus,input.seqs=Gekko.japonicus.genome_url,output.path="./Gekko.japonicus.tblastx.best.hits_seqs.fas")
Pogona.vitticeps.best.hits.seqs             <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Pogona.vitticeps,input.seqs=Pogona.vitticeps.genome_url,output.path="./Pogona.vitticeps.tblastx.best.hits_seqs.fas")
Crotalus.horridus.best.hits.seqs            <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.horridus,input.seqs=Crotalus.horridus.genome_url,output.path="./Crotalus.horridus.tblastx.best.hits_seqs.fas")
Crotalus.mitchellii.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Crotalus.mitchellii,input.seqs=Crotalus.mitchellii.genome_url,output.path="./Crotalus.mitchellii.tblastx.best.hits_seqs.fas")
Ophiophagus.hannah.best.hits.seqs           <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Ophiophagus.hannah,input.seqs=Ophiophagus.hannah.genome_url,output.path="./Ophiophagus.hannah.tblastx.best.hits_seqs.fas")
Pantherophis.guttatus.best.hits.seqs        <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Pantherophis.guttatus,input.seqs=Pantherophis.guttatus.genome_url,output.path="./Pantherophis.guttatus.tblastx.best.hits_seqs.fas")
Protobothrops.mucrosquamatus.best.hits.seqs <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Protobothrops.mucrosquamatus,input.seqs=Protobothrops.mucrosquamatus_genome.url,output.path="./Protobothrops.mucrosquamatus.tblastx.best.hits_seqs.fas")
Python.bivittatus.best.hits.seqs            <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Python.bivittatus,input.seqs=Python.bivittatus.genome_url,output.path="./Python.bivittatus.tblastx.best.hits_seqs.fas")
Viperus.berus.best.hits.seqs                <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Viperus.berus,input.seqs=Viperus.berus.genome_url,output.path="./Viperus.berus.tblastx.best.hits_seqs.fas")
Thamnophis.sirtalis.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Thamnophis.sirtalis,input.seqs=Thamnophis.sirtalis.genome_url,output.path="./Thamnophis.sirtalis.tblastx.best.hits_seqs.fas")
```

The output sequences were saved in fasta format and can be downloaded here: [Anolis.carolinensis.best.hits_seqs.fas](https://osf.io/efu2n/download), [Gekko.japonicus.best.hits_seqs.fas](https://osf.io/rpzvh/download), [Pogona.vitticeps.best.hits_seqs.fas](https://osf.io/hg9pq/download), [Crotalus.horridus.best.hits_seqs.fas](https://osf.io/9qbrx/download), [Crotalus.mitchellii.best.hits_seqs.fas](https://osf.io/q9ytz/download), [Ophiophagus.hannah.best.hits_seqs.fas](https://osf.io/ckxm9/download), [Pantherophis.guttatus.best.hits_seqs.fas](https://osf.io/9wd4a/download), [Protobothrops.mucrosquamatus.best.hits_seqs.fas](https://osf.io/hpdw9/download), [Python.bivittatus.best.hits_seqs.fas](https://osf.io/b3d6w/download), [Viperus.berus.best.hits_seqs.fas](https://osf.io/5tnr4/download), [Thamnophis.sirtalis.best.hits_seqs.fas](https://osf.io/psvg8/download).

6. I used the makeStatsTable function (REEs package) to perform multiple sequence alignment (MAFFT algorithm) for the best matches to each query exon sequence and to create a table of alignment statistics for each alignment. The information and statistics in the output table are summarized in Table 3.

```
### Create a character vector holding the paths to each of the ".*.best.hits_seqs.fas" files generated in step 5. 
input.seqs.paths <- c("Anolis.carolinensis.best.hits_seqs.fas", "Gekko.japonicus.best.hits_seqs.fas", "Pogona.vitticeps.best.hits_seqs.fas", "Crotalus.horridus.best.hits_seqs.fas", "Crotalus.mitchellii.best.hits_seqs.fas", "Ophiophagus.hannah.best.hits_seqs.fas", "Pantherophis.guttatus.best.hits_seqs.fas", "Protobothrops.mucrosquamatus.best.hits_seqs.fas", "Python.bivittatus.best.hits_seqs.fas", "Viperus.berus.best.hits_seqs.fas", "Thamnophis.sirtalis.best.hits_seqs.fas")

### Align homologous sequences and make a table to summarize data in each alignment (one row per aligned locus).
stats.table.all  <- makeStatsTable(input.seqs=input.seqs.paths, input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp, output.path="./statsTable_REEs_SnakeCap.txt", species=c("Anolis carolinensis","Gekko japonicus","Pogona vitticeps","Crotalus horridus","Crotalus mitchellii","Ophiophagus hannah","Pantherophis guttatus", "Protobothrops mucrosquamatus", "Python bivittatus","Vipera berus","Thamnophis sirtalis"), subgroup=c("Crotalus horridus","Crotalus mitchellii","Ophiophagus hannah","Pantherophis guttatus", "Protobothrops mucrosquamatus", "Python bivittatus","Vipera berus","Thamnophis sirtalis"),alignments.out=alignments.dir, species.gff=11)
```
The output table [statsTable_REEs_SnakeCap.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/statsTable_REEs_SnakeCap.txt) includes summary statistics and information for 66,489 alignments (**XXXXX.tar.xz**).

**Table 3.** Description of information and statistics in each column of the output table generated by the function makeStatsTable in step 6. If you are using this table to interpret your own results with a different number of input species, then the column number will differ from what is shown here. Also, column names that include a species name will differ if you used different input species. The species name in columns 1 and 7 will always be the name of the primary species chosen in step 1.
column number|column name|column description
---|---|---
1|Thamnophis_sirtalis.locus|Contig accession and start and end position of the query sequence.
2|num.Species|Number of species (individuals) in the alignment
3|CountCover|Number of sites with at least four individuals with non-missing or ambiguous data
4|absolutePIS|Number of parsimony informative sites
5|percentPIS|percent of the sites with at least four individuals with non-missing or ambiguous data that are parsimony informative (i.e., percent of sites that could be parsimony informative that actually are)
6|mean.pident|Mean pairwise percent of sites identical to *T. sirtalis*; this is the mean of columns 8–17.
7|pident.Thamnophis_sirtalis|Percent of sites identical between *T. sirtalis* query and *T. sirtalis* best match; this should always be 100 and any other value indicates the tblastx match that was chosen as the best match (step 4 reportBestMatches function) is not the true best match.
8|pident.Anolis_carolinensis|Percent of sites identical between *T. sirtalis* query and *A. carolinensis* best match.
9|pident.Gekko_japonicus|Percent of sites identical between *T. sirtalis* query and *G. japonicus* best match.
10|pident.Pogona_vitticeps|Percent of sites identical between *T. sirtalis* query and *P. vitticeps* best match.
11|pident.Crotalus_horridus|Percent of sites identical between *T. sirtalis* query and *C. horridus* best match.
12|pident.Crotalus_mitchellii|Percent of sites identical between *T. sirtalis* query and *C. mitchellii* best match.
13|pident.Ophiophagus_hannah|Percent of sites identical between *T. sirtalis* query and *O. hannah* best match.
14|pident.Pantherophis_guttatus|Percent of sites identical between *T. sirtalis* query and *P. guttatus* best match.
15|pident.Protobothrops_mucrosquamatus|Percent of sites identical between *T. sirtalis* query and *P. mucrosquamatus* best match.
16|pident.Python_bivittatus|Percent of sites identical between *T. sirtalis* query and *P. bivittatus* best match.
17|pident.Viperus_berus|Percent of sites identical between *T. sirtalis* query and *V. berus* best match.
18|gene.name|Name of the gene of the *T. sirtalis* query locus. The is the name associated with the query sequence's CDS feature in the *T. sirtalis* GFF genome feature table.
19|locus.length.Thamnophis_sirtalis|Length of the query sequence
20|mean.variable.sites| Mean pairwise fraction of sites variable compared to *T. sirtalis*. This is simply (100-mean.pident)/100
21|min.pident.all|Minimum of percent of sites identical to *T. sirtalis* among the species. In other words, the minimum of columns 8–17. This would be calculated across different columns if different species were used.
22|min.pident.subgroup|Minimum of percent of sites identical to *T. sirtalis* among the snake species compared. In other words, the minimum of columns 11–17. The would be calculated across different columns if the subgroup was defined differently.
23|alignment.width|Width of the alignment. This column was not generated for the SnakeCap dataset because I used an older version of makeStatsTable, which did not include this column in the output table.

  7. I used the function pick.loci (REEs package) to choose a set of REEs for sequence capture. I applied the following criteria to select REEs: (1) I filtered out loci if minimum percent genetic similarity (mean among snakes) to *T. sirtalis* was < 65% or = 100%; (2) for genes with multiple exons, I kept only the exon with the lowest mean pairwise genetic distance to *T. sirtalis* (16,650 exons pass this step); (3) the maximum number of baits required to target the exons was set to 20,000 (this constraint was imposed by the 20K my-baits kit) with bait size = 120nt and 50% bait tiling; (4) the maximum number of nucleotides to target was set to 1.2Mb (also constrained by the my-baits kit); (5) of the sets of exons that meet criteria 1–4, I chose the set that maximizes the total number of variable sites relative to *T. sirtalis*.

<!--
**The number of nucleotides targetted should be determined by the number of baits in the kit, length of baits, and tiling amount...For each row in the input stats matrix, calculate how many baits would be needed to target the exon using the bait length and tiling amount, and then use rollSum function on the number of baits (after other filtering steps and sorting steps are completed). Keep rows in which rollSum(number.of.baits) < 20,000. However, check if the 1.2Mb threshold wasn't determined by Arbor...
-->

```
### Read the table generated by makeStatsTable function
stats.table.all <- data.table::fread("./statsTable_REEs_SnakeCap.txt",header=T)

### Filter the stats table to the optimal set of REEs
stats.table.best <- REEs::pick.loci(statsTable.path=stats.table.all,primary.species="Thamnophis sirtalis", output.path="./stats_data_FastestExonPerGene_best_28Dec2020.tsv", species.subgroup=c("Crotalus horridus","Crotalus mitchellii","Ophiophagus hannah","Pantherophis guttatus", "Protobothrops mucrosquamatus", "Python bivittatus","Vipera berus","Thamnophis sirtalis"),pident.keep=c(65,100),max.loci.per.gene=1, min.num.species="all",max.capture.coverage=1200500,fast.stat="pident",use.min.pident.subgroup=T)
```

<!---
### Note 1: Using the code above retains 2,068 REEs in stats.table.best, of which 2,067 are the same as those previously picked and included in the file "stats_data_FastestExonPerGene_best.tsv". The output table stats.table.best includes "NW_013658076.1:768357-769730", which was not previously selected for inclusion in "stats_data_FastestExonPerGene_best.tsv"; conversely, "stats_data_FastestExonPerGene_best.tsv" includes "NW_013658076.1:768360-769730", "NW_013657914.1:650880-651587", "NW_013662230.1:5224-5446", and "NW_013659343.1:156011-156388", which are not included in the output table stats.table.best
### Note 2: "NW_013658076.1:768357-769730" of stats.table.best is essentially the same locus as "NW_013658076.1:768360-769730" in stats_data_FastestExonPerGene_best.tsv; The difference is because the latest version of the REEs::reportBestMatches function filters matches that are subsequences of other matches.
### Note 3: running pick.loci with max.capture.coverage=1200500 (rather than 1200000) includes the loci "NW_013657914.1:650880-651587" and "NW_013662230.1:5224-5446", which were included in "stats_data_FastestExonPerGene_best.tsv", but "NW_013659343.1:156011-156388" is still missing.
--->

The output table includes 2,070 REEs and can be downloaded here [stats_data_FastestExonPerGene_best_28Dec2020.tsv](https://github.com/JeffWeinell/SnakeCap/raw/main/REEs/stats_data_FastestExonPerGene_best_28Dec2020.tsv). The format of the output table is the same as the format of the input table; see Table 3 for column descriptions.

<!--
An updated version of this stats table stats_data_FastestExonPerGene_best.tsv that includes the WeinellEntry locus names is [stats_data_FastestExonPerGene_best_20Nov2020.tsv](https://github.com/JeffWeinell/SnakeCap/raw/main/REEs/stats_data_FastestExonPerGene_best_20Nov2020.tsv).
-->

8. Expand the target region to include noncoding DNA upstream and downstream of the loci identified in step 7. The amount of noncoding DNA included in each expanded target is a function of the bait length (120nt) and exon length, such that the expanded target length is a multiple of the bait length.

```
### Calculate start and end coordinates of expanded targets (REEs + flanking noncoding regions) such that the expanded target length is a multiple of the bait length.

### Read the output table of the function pick.loci in step 7.
stats.table.best       <- as.data.frame(data.table::fread("stats_data_FastestExonPerGene_best_28Dec2020.tsv",header=T,sep="\t"))

### Define a vector holding the locus lengths
locus.lengths          <- stats.table.best[,"locus.length.Thamnophis_sirtalis"]

### Subtract one from each value in locus.lengths
### NOTE: this is done for reproducibility of SnakeCap methods,
### but skip the next line when designing baits for your own project.
locus.lengths <- locus.lengths-1

### Create a numeric matrix of the start and end coordinates of the REEs. This information is extracted from the first column of the input table.
REEs.coordinates       <- mat.strsplit(gsub(".*:","",stats.table.best[,"Thamnophis_sirtalis.locus"]),split="-")
mode(REEs.coordinates) <- "numeric"

### Define the length of baits
bait.length    <- 120

### Calculate lengths for the expanded targets (= locus.lengths rounded up to the nearest multiple of bait.length)
target.lengths <- ceiling(locus.lengths/bait.length)*bait.length

### Calculate amount of upstream noncoding DNA to include
upstream.lengths   <- floor(((target.lengths-locus.lengths)/2))

### Calculate amount of downstream noncoding DNA to include
downstream.lengths <- target.lengths-(locus.lengths+upstream.lengths)

### Define coordinates for the expanded targets (contig accession ID, start position, end position). The contig accession IDs are extracted from the first column of the input table. 
targets.contig <- mat.strsplit(gsub(":.*","",stats.table.best[,"Thamnophis_sirtalis.locus"]),split="-")
targets.start  <- REEs.coordinates[,1]-upstream.lengths
targets.end    <- REEs.coordinates[,2]+downstream.lengths

### Define URL to the T. sirtalis genome sequences. This URL is stored as a character string in the REEs package, and can be accessed with the datasets() function.
Thamnophis.sirtalis_genome.url <- REEs::datasets(1)[1,2]

### Extract and save sequences for the expanded targets
### When selecting targets, I suggest setting the trim.ambiguous parameter to TRUE when using the get_ncbi_sequences function (REEs package); this parameter was not implemented in the SnakeCap project, and, therefore, for the sake of reproducibility, trim.ambiguous is set to FALSE here.

REEs.expanded <- get_ncbi_sequences(outfile="./REEs.expanded.fas",input.seqs=Thamnophis.sirtalis_genome.url, accessionList=targets.contig, startList= targets.start, endList=targets.end, if.outside.range="partial",trim.ambiguous = FALSE)
```

The output sequences (expanded REEs targets) can be downloaded here: [REEs.expanded.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/REEs/REEs.expanded.fas).

Results:
Nine of the expanded targets were only partially expanded, meaning that their sequence lengths (actually sequence lengths minus one, because of a bug in the code that has since been fixed) were not a multiple of the bait length. Three of these targets were located near the end of the reference contig and another six had terminal strings of ambiguous bases (Ns) that were trimmed (Table 4). Note: 57 other targets that were fully expanded had terminal Ns, but these Ns were not trimmed because they were not detected. The latest version of the get_ncbi_sequences function (REEs package) has an option to trim terminal Ns, but this option was not present when SnakeCap targets were chosen.

Additionally, NW_013658076.1:768325-769765 (WeinellEntry1658) was targetted instead of NW_013658076.1:768324-769764; these are nearly identical (shifted by only one base on the contig) and I am not sure why this change was made.

Five duplicate pairs of REEs (each pair with identical sequences) were present in the output of step 8 (Table 5). Only one of these pairs was recognized/identified (and filtered manually) prior to submitting target sequences to Arbor Biosciences.

The remaining 2,068 REEs (expanded targets) were submitted to Arbor Biosciences for ulstrastringent filtering and probe design, and can be downloaded here: [REEs.expanded.final.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/REEs/REEs.expanded.final.fas). These were actually submitted to Arbor in two batches: [Version1_Target-loci_Jeff-Weinell_10Sep2018.fasta](https://github.com/JeffWeinell/SnakeCap/raw/main/ArborFiles/Version1_Target-Loci_Jeff-Weinell_10Sep2018.fasta) and [Version2_additional-targets_20Sep2018.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/ArborFiles/Version2_additional-targets_Entry1899to3152_20Sep2018.fasta). See [ultra-stringent filtering](#ultrastringentFiltering) section.

Arbor performed ultrastringent filtration on the 2,068 REEs retained from step 8 (after removing two identical sequences; pair 1 Table 5). Ultrastringent filtering resulted in the removal of 203 REEs (1,865 REEs retained). Of the 203 REEs that were filtered, 76 were filtered because no baits could be designed for these loci (70 of these are listed in [Version1-loci-removed_ZeroBaitCoverageLoci.tsv](https://git.io/JLiEu) and six are listed in [Version2-loci-removed_ZeroBaitCoverageLoci.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version2-loci-removed_ZeroBaitCoverageLoci.tsv)); and 127 REEs were filtered because all proposed baits were non-specific within the *T. sirtalis* genome (eight of these are listed in [Version1-loci-removed_baits-nonspecific.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version1-loci-removed_nonspecific-baits.tsv) and 119 are listed in [Version2-loci-removed_baits-nonspecific.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version2-loci-removed_baits-nonspecific.tsv)).

Of the 1,865 REEs that passed ultrastringent filtering, 212 were removed to allow a fraction of the 20K baits to be used to target other types of loci (UCEs, MHC genes, scalation genes, vision genes, and ddRAD-like loci), and these removed REEs are listed in the file [Version3-loci-removed_others.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version3-loci-removed_others.tsv). Baits for the remaining 1,653 REEs were synthesized by Arbor (mybaits 20K bait kit: product no. 3001160).

Table 4. The nine partially expanded targets from output of step 8 and reason why targets were partially rather than fully expanded.
WeinellEntry name|Contig accession ID|Nucleotide range of partially expanded target|Nucleotide range that would have been targetted if it had been possible|Reason why target sequence partially rather than fully expanded.
---|---|---|---|---
WeinellEntry959 |NW_013657802.1|1404174-1405963 |1404174-1405974|contig length 140,5963 nt
WeinellEntry1800|NW_013662380.1|5926-6377       |5926-6406      |contig length 6,377 nt
WeinellEntry2040|NW_013661694.1|23750-23983     |23750-23990    |contig length 23,983 nt
WeinellEntry2147|NW_013657804.1|885274-886205   |885274-886234  |886206-886234 all Ns
WeinellEntry2148|NW_013659850.1|163196-163554   |163196-163556  |163555-163556 all Ns
WeinellEntry2149|NW_013659217.1|238701-238933   |238701-238941  |238933-238941 all Ns
WeinellEntry2151|NW_013658527.1|12230-12575     |12230-12590    |12576-12590 all Ns
WeinellEntry2150|NW_013658733.1|434060-434488   |434008-434488  |434008-434059 all Ns
WeinellEntry2152|NW_013658527.1|14050-14395     |14035-14395    |14035-14049 all Ns

Table 5. Pairs of REEs having identical sequences that were included in the ouput of the pick.loci function (step 8). These were subsequently filtered, either immediately before or after application of Arbor's ultrastringent filtering algorithm. The latest version of the pick.loci function has an option to filter REEs if the bitscores of the top matches are too similar according to a user-defined threshold.
Contig Accession ID|Start Position|End Position|Sequence/Pair ID|Other ID|Step when filtered
---|---|---|---|---|---
NW_013657725.1|467272|467752|1||manually, after using pick.loci function and before ultrastringent filtering
NW_013657725.1|516435|516915|1||manually, after using pick.loci function and before ultrastringent filtering
NW_013659646.1|15015|15975|2|WeinellEntry44|during ultrastringent filtering
NW_013659646.1|24066|25026|2|WeinellEntry45|during ultrastringent filtering
NW_013657804.1|817033|817753|3|WeinellEntry496|during ultrastringent filtering
NW_013657804.1|820482|821202|3|WeinellEntry497|during ultrastringent filtering
NW_013658610.1|34634|34874|4|WeinellEntry589|during ultrastringent filtering
NW_013658610.1|59790|60030|4|WeinellEntry590|during ultrastringent filtering
NW_013658165.1|745869|746109|5|WeinellEntry1670|during ultrastringent filtering
NW_013658165.1|748461|748701|5|WeinellEntry1671|during ultrastringent filtering

<!--
REEs that failed ultrastringent filtration: 70
REEs removed because all baits were non-specific: 123
REEs removed to make room for other loci: 212
Total REEs removed: 70+123+212 = 405
1,653 REEs were synthesized and included in the bait kit.
10 REEs unaccounted for...
-->

<!--- Next step was done but wasn't necessary. It might be good to run this as a sanity check when picking REEs for future probe sets
7. I used the R function align.and.concatenate.best.exons to perform multiple sequence alignment (MAFFT algorithm) for each exon in the "stats_data_FastestExonPerGene_best.txt" file (which was generated by the pick.loci function in step 5), and to concatenate exon alignments and generate an associated partition file. Single-locus alignments, the concatenated loci alignment, and the partition file were each saved to file. I estimated gene trees from these alignments (using IQTREE) to get a sense for how much phylogenetic information each locus contained.
--->

<!--- Idk what I need these commented out lines for:
8. For each locus, I chose the target region that will be used for probe design by considering the length of each target exon, probe length (120bp), and tiling regime (50% overlap between adjacent probes).
Li = length of ith target exon
nprobes,i = minimum # of probes needed to cover Li
Lprobes,i = # of nucleotides covered by nprobes,i
nflank,i = # of exon-flanking nucleotides targeted (total, upstream + downstream)
nflank.5',i = # of nucleotides targeted upstream of exon
= (nflank,i)/2 rounded down to nearest integer
nflank.3',i = # of nucleotides targeted downstream of exon
= nflank,i - nflank.5',i
--->

<a name="Methods.SelectingUCEs"></a>
## Selecting the set of target UCEs

<!-- <a name="Methods.SelectingUCEs.overview"></a> -->
#### Overview:

Target UCEs include 907 of the 3,260 UCEs previously identified in *Micrurus fulvius* (Streicher and Wiens, 2017; **Table X**). First, I filtered the full set of *Micrurus* UCEs to only include those present in all NCBI snake genomes (n = 2,968 UCEs). Then, I filtered the shared set of UCEs to only include those with an alignment width > 200nt (2,551 UCEs retained; each UCE alignment included the eight snakes with published genomes). Next, I removed the following UCEs: uce-1843, uce-2179, uce-2433, uce-2465, uce-2498, uce-2890, uce-2960, and uce-3354 (not sure why I did this yet). I sorted the remaining 2,543 UCEs by UCE name and retained the first 1,000 loci in this set. Arbor Biosciences was able to synthesize probes for 907 of the 1,000 proposed target UCEs.

<!-- <a name="Methods.SelectingUCEs.detailed"></a> -->
#### Details:

1. I downloaded the set of *Micurus fulvius* UCEs (n = 3,260) identified by Streicher and Wiens (2017). These were available as a fasta file called [micrurus_UCEs.fa](https://git.io/JLiu2).

2. To find *Micrurus* UCEs in each of the other snake genomes, I queried each *Micrurus* UCE against each snake genome using blastn algorithm (saving ≤ 50 matches per query), which was implemented using the R wrapper function REEs::blast.

```
# Runs BLASTN
Crotalus.horridus.UCEs.50hits            <- REEs::blast(method="blastn",subject=Crotalus.horridus.genome_url, query="./micrurus_UCEs.fa",table.out="./Crotalus.horridus.blastn.UCEs.50hits.txt")
Crotalus.mitchellii.UCEs.50hits          <- REEs::blast(method="blastn",subject=Crotalus.mitchellii.genome_url, query="./micrurus_UCEs.fa",table.out="./Crotalus.mitchellii.blastn.UCEs.50hits.txt")
Ophiophagus.hannah.UCEs.50hits           <- REEs::blast(method="blastn",subject=Ophiophagus.hannah.genome_url, query="./micrurus_UCEs.fa",table.out="./Ophiophagus.hannah.blastn.UCEs.50hits.txt")
Pantherophis.guttatus.UCEs.50hits        <- REEs::blast(method="blastn",subject=Pantherophis.guttatus.genome_url, query=Thamnophis.sirtalis_exome,table.out="./Pantherophis.guttatus.blastn.UCEs.50hits.txt")
Protobothrops.mucrosquamatus.UCEs.50hits <- REEs::blast(method="blastn",subject=Protobothrops.mucrosquamatus.genome_url, query="./micrurus_UCEs.fa",table.out="./Protobothrops.mucrosquamatus.blastn.UCEs.50hits.txt")
Python.bivittatus.UCEs.50hits            <- REEs::blast(method="blastn",subject=Python.bivittatus.genome_url, query="./micrurus_UCEs.fa",table.out="./Python.bivittatus.blastn.UCEs.50hits.txt")
Viperus.berus.UCEs.50hits                <- REEs::blast(method="blastn",subject=Viperus.berus.genome_url, query="./micrurus_UCEs.fa",table.out="./Viperus.berus.blastn.UCEs.50hits.txt")
Thamnophis.sirtalis.UCEs.50hits          <- REEs::blast(method="blastn",subject=Thamnophis.sirtalis.genome_url, query="./micrurus_UCEs.fa",table.out="./Thamnophis.sirtalis.blastn.UCEs.50hits.txt")
```

Output tables from BLASTN: **Crotalus.horridus.blastn.UCEs.50hits.txt**, **Crotalus.mitchellii.blastn.UCEs.50hits.txt**, **Ophiophagus.hannah.blastn.UCEs.50hits.txt**, **Pantherophis.guttatus.blastn.UCEs.50hits.txt**, **Protobothrops.mucrosquamatus.blastn.UCEs.50hits.txt**, **Python.bivittatus.blastn.UCEs.50hits.txt**, **Thamnophis.sirtalis.blastn.UCEs.50hits.txt**, and **Vipera.berus.blastn.UCEs.50hits.txt**.

 3. I used the function reportBestMatches to filter the hit tables to include only the best match for each query UCE.

```
UCEs.best.hits.Crotalus.horridus             <- REEs::reportBestMatches(input.table=Crotalus.horridus.blastn.UCEs.50hits, output.table.path="Crotalus.horridus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Crotalus.mitchellii           <- REEs::reportBestMatches(input.table=Crotalus.mitchellii.blastn.UCEs.50hits, output.table.path="Crotalus.mitchellii.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Ophiophagus.hannah            <- REEs::reportBestMatches(input.table=Ophiophagus.hannah.blastn.UCEs.50hits, output.table.path="Ophiophagus.hannah.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Pantherophis.guttatus         <- REEs::reportBestMatches(input.table=Pantherophis.guttatus.blastn.UCEs.50hits, output.table.path="Pantherophis.guttatus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Protobothrops.mucrosquamatus  <- REEs::reportBestMatches(input.table=Protobothrops.mucrosquamatus.blastn.UCEs.50hits, output.table.path="Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Python.bivittatus             <- REEs::reportBestMatches(input.table=Python.bivittatus.blastn.UCEs.50hits, output.table.path="Python.bivittatus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Thamnophis.sirtalis           <- REEs::reportBestMatches(input.table=Thamnophis.sirtalis.blastn.UCEs.50hits, output.table.path="Thamnophis.sirtalis.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Viperus.berus                 <- REEs::reportBestMatches(input.table=Vipera.berus.blastn.UCEs.50hits, output.table.path="Vipera.berus.blastn.UCEs.best.hits.txt")
```

Output tables were saved to the files: [Crotalus.horridus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.horridus.blastn.UCEs.best.hits.txt), [Crotalus.mitchellii.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.mitchellii.blastn.UCEs.best.hits.txt), [Ophiophagus.hannah.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus.hannah.blastn.UCEs.best.hits.txt) , [Pantherophis.guttatus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis.guttatus.blastn.UCEs.best.hits.txt), [Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt), [Python.bivittatus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Python.bivittatus.blastn.UCEs.best.hits.txt), [Thamnophis.sirtalis.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis.sirtalis.blastn.UCEs.best.hits.txt), [Vipera.berus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Vipera.berus.blastn.UCEs.best.hits.txt).

4. I used the function REEs::get.seqs.from.blastTable to extract and save the set of best-match UCEs from each genome.
```
#### Read in the tables of best matches generated from step 3 (if these are not already loaded).
UCEs.best.hits.Crotalus.horridus             <- as.data.frame(data.table::fread("./Crotalus.horridus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Crotalus.mitchellii           <- as.data.frame(data.table::fread("./Crotalus.mitchellii.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Ophiophagus.hannah            <- as.data.frame(data.table::fread("./Ophiophagus.hannah.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Pantherophis.guttatus         <- as.data.frame(data.table::fread("./Pantherophis.guttatus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Protobothrops.mucrosquamatus  <- as.data.frame(data.table::fread("./Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Python.bivittatus             <- as.data.frame(data.table::fread("./Python.bivittatus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Thamnophis.sirtalis           <- as.data.frame(data.table::fread("./Thamnophis.sirtalis.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Viperus.berus                 <- as.data.frame(data.table::fread("./Vipera.berus.blastn.UCEs.best.hits.txt"),header=T)

#### Define the genome URL paths. The URLs for genomes used in the SnakeCap study can be called using the datasets function of REEs package.
Crotalus.horridus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus horridus"),2]
Crotalus.mitchellii.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus mitchellii"),2]
Ophiophagus.hannah.genome_url           <- REEs::datasets(1)[which(datasets(1)[,1]=="Ophiophagus hannah"),2]
Pantherophis.guttatus.genome_url        <- REEs::datasets(1)[which(datasets(1)[,1]=="Pantherophis guttatus"),2]
Protobothrops.mucrosquamatus.genome_url <- REEs::datasets(1)[which(datasets(1)[,1]=="Protobothrops mucrosquamatus"),2]
Python.bivittatus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Python bivittatus"),2]
Viperus.berus.genome_url                <- REEs::datasets(1)[which(datasets(1)[,1]=="Viperus berus"),2]
Thamnophis.sirtalis.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Thamnophis sirtalis"),2]

#### Extracts the sequence of the best match of each UCE from each snake genome.
Crotalus.horridus.best.hits.seqs            <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Crotalus.horridus, input.seqs=Crotalus.horridus.genome_url, output.path="./Crotalus_horridus_UCEs.fasta")
Crotalus.mitchellii.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Crotalus.mitchellii, input.seqs=Crotalus.mitchellii.genome_url, output.path="./Crotalus_mitchellii_UCEs.fasta")
Ophiophagus.hannah.best.hits.seqs           <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Ophiophagus.hannah, input.seqs=Ophiophagus.hannah.genome_url, output.path="./Ophiophagus_hannah_UCEs.fasta")
Pantherophis.guttatus.best.hits.seqs        <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Pantherophis.guttatus, input.seqs=Pantherophis.guttatus.genome_url, output.path="./Pantherophis_guttatus_UCEs.fasta")
Protobothrops.mucrosquamatus.best.hits.seqs <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Protobothrops.mucrosquamatus, input.seqs=Protobothrops.mucrosquamatus_genome.url, output.path="./Protobothrops_mucrosquamatus_UCEs.fasta")
Python.bivittatus.best.hits.seqs            <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Python.bivittatus, input.seqs=Python.bivittatus.genome_url, output.path="./Python_bivittatus_UCEs.fasta")
Viperus.berus.best.hits.seqs                <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Viperus.berus, input.seqs=Viperus.berus.genome_url, output.path="./Viperus_berus_UCEs.fasta")
Thamnophis.sirtalis.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Thamnophis.sirtalis, input.seqs=Thamnophis.sirtalis.genome_url, output.path="./Thamnophis_sirtalis_UCEs.fas")
```

Output sequences in fasta format can be downloaded here: [Crotalus.horridus.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.horridus.UCEs.fas), [Crotalus.mitchellii.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.mitchellii.UCEs.fas), [Ophiophagus.hannah.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus.hannah.UCEs.fas), [Pantherophis.guttatus.UCEs.fasta](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis.guttatus.UCEs.fas), [Protobothrops.mucrosquamatus.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops.mucrosquamatus.UCEs.fas), [Python.bivittatus.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Python.bivittatus.UCEs.fas), [Thamnophis.sirtalis.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis.sirtalis.UCEs.fas), [Vipera.berus.UCEs.fasta](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Vipera.berus.UCEs.fas)

Important note: Sequence headers in output fasta files (i.e., the text on lines starting with ">") have the format "querySequenceName_Subject=subjectContigName_ContigStart_ContigEnd". For example, the *Crotalus horridus* sequence that is the best match to UCE 1003 of *Micrurus fulvius* has the header "micrurus_fulvius_uce-1003_Subject=LVCR01005387.1_10818_10626". This can be a bit confusing at first because the only species mentioned in the sequence header is Micrurus fulvius, even though the sequence is actually from *Crotalus horridus*.

To avoid confusion, I saved a copy of each species' UCE sequences after renaming sequences with the more intuitive header format "Genus_species_uce-XXXX", where XXXX is a number. For example: "Crotalus_horridus_uce-1003" is used as the name for the homolog of "micrurus_fulvius_uce-1003".

```
### Renaming each species' UCE sequences to match the format of sequence names used by Streicher and Wiens (2017): "Genus_species_uce-XXXX", where XXXX is a number.

names(Crotalus.horridus.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Crotalus_horridus_",""),names(Crotalus.horridus.best.hits.seqs))
names(Crotalus.mitchellii.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Crotalus_mitchellii_",""),names(Crotalus.mitchellii.best.hits.seqs))
names(Ophiophagus.hannah.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Ophiophagus_hannah_",""),names(Ophiophagus.hannah.best.hits.seqs))
names(Pantherophis.guttatus.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Pantherophis_guttatus_",""),names(Pantherophis.guttatus.best.hits.seqs))
names(Protobothrops.mucrosquamatus.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Protobothrops_mucrosquamatus_",""),names(Protobothrops.mucrosquamatus.best.hits.seqs))
names(Python.bivittatus.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Python_bivittatus_",""),names(Python.bivittatus.best.hits.seqs))
names(Viperus.berus.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Viperus_berus_",""),names(Viperus.berus.best.hits.seqs))
names(Thamnophis.sirtalis.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Thamnophis_sirtalis_",""),names(Thamnophis.sirtalis.best.hits.seqs))

### Save a copy of each of each species' renamed UCE sequences
writeXStringSet(Crotalus.horridus.best.hits.seqs,"Crotalus.horridus.UCEs_renamed.fas")
writeXStringSet(Crotalus.mitchellii.best.hits.seqs,"Crotalus.mitchellii.UCEs_renamed.fas")
writeXStringSet(Ophiophagus.hannah.best.hits.seqs,"Ophiophagus.hannah.UCEs_renamed.fas")
writeXStringSet(Pantherophis.guttatus.best.hits.seqs,"Pantherophis.guttatus.UCEs_renamed.fas")
writeXStringSet(Protobothrops.mucrosquamatus.best.hits.seqs,"Protobothrops.mucrosquamatus.UCEs_renamed.fas")
writeXStringSet(Python.bivittatus.best.hits.seqs,"Python.bivittatus.UCEs_renamed.fas")
writeXStringSet(Viperus.berus.best.hits.seqs,"Thamnophis.sirtalis.UCEs_renamed.fas")
writeXStringSet(Thamnophis.sirtalis.best.hits.seqs,"Vipera.berus.UCEs_renamed.fas")
```

The renamed sequences can be downloaded here: [Crotalus.horridus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.horridus.UCEs_renamed.fas), [Crotalus.mitchellii.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.mitchellii.UCEs_renamed.fas), [Ophiophagus.hannah.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus.hannah.UCEs_renamed.fas), [Pantherophis.guttatus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis.guttatus.UCEs_renamed.fas), [Protobothrops.mucrosquamatus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops.mucrosquamatus.UCEs_renamed.fas), [Python.bivittatus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Python.bivittatus.UCEs_renamed.fas), [Thamnophis.sirtalis.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis.sirtalis.UCEs_renamed.fas), [Vipera.berus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Vipera.berus.UCEs_renamed.fas)

<!---
5. I used the function REEs::align.bestHit.UCEs to identify and align the set of UCEs found in all snake genomes. This function invokes MAFFT to perform multisequence alignment.
```
align.bestHit.UCEs(species.UCEs.filepaths=list.files(path="~/UCEs.In.Snake.Genomes/",full.names=T), output.dir="~/MAFFT-aligned-UCEs", species=c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus", "Protobothrops_mucrosquamatus","Pantherophis_guttatus"))
```
--->

5. I used the makeStatsTable function (REEs package) to perform multiple sequence alignment (MAFFT algorithm) for the best matches to each UCE, and to create a table holding some statistics for each alignment.

```
### Create a character vector holding the paths to the fasta files generated in step 5. 
input.seqs.UCEs.paths <- c("Crotalus_horridus_UCEs.fasta", "Crotalus_mitchellii_UCEs.fasta", "Ophiophagus_hannah_UCEs.fasta", "Pantherophis_guttatus_UCEs.fasta", "Protobothrops_mucrosquamatus_UCEs.fasta", "Python_bivittatus_UCEs.fasta", "Vipera_berus_UCEs.fasta", "Thamnophis_sirtalis_UCEs.fas")

### Define output directory for UCE alignments
UCEs.alignments.dir <- "~/MAFFT-aligned-UCEs"

### Align sequences for each UCE and make a table to summarize data in each alignment (one row per UCE locus).
stats.table.UCEs.all  <- makeStatsTable(input.seqs=input.seqs.UCEs.paths, input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp, output.path="./statsTable_REEs_SnakeCap.txt", species=c("Crotalus horridus","Crotalus mitchellii","Ophiophagus hannah","Pantherophis guttatus", "Protobothrops mucrosquamatus", "Python bivittatus","Vipera berus","Thamnophis sirtalis"),alignments.out=UCEs.alignments.dir, species.gff=8)
```

6. Final size filtering, sorting, and UCE selection steps were performed in R. *Thamnophis sirtalis* sequences for the selected UCEs (n = 1,000) were submitted to Arbor Biosciences for probe design. Requires Biostrings and ape packages.

```
UCE.alignment.filenames    <- list.files(path="~/MAFFT-aligned-UCEs",full.names=T)
UCE.shortnames             <- mgsub(patt=c(".fasta","uce-"),repl=c("","UCE."),subj=list.files(path="~/MAFFT-aligned-UCEs",full.names=F))

# Reads in UCE alignments.
for(i in 1:length(UCE.shortnames)){                                                           
	assign(x=UCE.shortnames[i],value=readDNAStringSet(filepath=UCE.alignment.filenames[i]))
}

# Hold UCE alignments in a list, sorted by UCE name
alignments.sorted <- mget(sort(UCE.shortnames))

### Renaming the individuals in each alignment (because names were truncated when saving alignments in the previous step).
for(i in 1:length(alignments.sorted)){
	names(alignments.sorted[[i]]) <- c("Thamnophis_sirtalis", "Crotalus_horridus", "Protobothrops_mucrosquamatus", "Ophiophagus_hannah", "Vipera_berus", "Crotalus_mitchellii", "Pantherophis_guttatus", "Python_bivittatus")
}

### Calculate mean pairwise genetic distance among individuals in each UCE alignment
mean.pdist <- vector(mode="numeric",length=length(alignments))
for(i in 1:length(alignments.sorted)){
	mean.pdist[i] <- round(mean(dist.dna(as.DNAbin(alignments.sorted[[i]]),model="raw",pairwise.deletion=T)),digits=3)
}
names(mean.pdist) <- names(alignments.sorted)
alignment.lengths <- unlist(lapply(alignments.sorted,FUN=function(X){unique(width(X))}))

# Remove UCE alignments if alignment width ≤ 200nt. Result: 417 UCEs removed and 2,551 retained
alignments.sorted.filtered <- alignments.sorted[which(alignment.lengths > 200)]

# Manually filter out "UCE.1843" "UCE.2179" "UCE.2433" "UCE.2465" "UCE.2498" "UCE.2890" "UCE.2960" "UCE.3354" "UCE.3501".
alignments.sorted.filtered[c("UCE.1843","UCE.2179","UCE.2433","UCE.2465","UCE.2498","UCE.2890","UCE.2960","UCE.3354","UCE.3501")] <- NULL

# Keep first 1,000 UCE alignments in alignments.sorted.filtered
alignments.sorted.filtered.1000 <- alignments.sorted.filtered[1:1000]

# Extract and save the Thamnophis sirtalis sequence (gaps removed) from each alignment in alignments.sorted.filtered.1000
library(DECIPHER) # need this package for the RemoveGaps function
target.UCEs <- list(); length(target.UCEs) <- length(alignments.sorted.filtered.1000)
for(i in 1:length(alignments.sorted.filtered.1000)){
	target.UCEs[[i]] <- alignments.sorted.filtered.1000[[i]]$Thamnophis_sirtalis
}
target.UCEs        <- RemoveGaps(DNAStringSet(target.UCEs))
names(target.UCEs) <- names(alignments.sorted.filtered.1000)
writeXStringSet(x=target.UCEs,filepath=<outputFilepath>,format="fasta")
```

7. Probes were designed for 907 of the 1,000 UCEs submitted to Arbor Biosciences.

<a name="Methods.SelectingddRAD"></a>
## Selecting the set of target ddRAD-like loci

#### Overview:

#### Details:

1. grep Sbfi recognition site in *T. baileyi* genome (sense strand contigs); output = three column hit table containing the "contig accession", "start position", "end position"
2. grep EcoRI recognition site in *T. baileyi* genome (sense strand contigs); output = three column hit table containing the "contig accession", "start position", "end position"
3. grep Sbfi recognition site in *T. baileyi* genome (antisense strand contigs); output = three column hit table containing the "contig accession", "start position", "end position"
4. grep EcoRI recognition site in *T. baileyi* genome (antisense strand contigs); output = three column hit table containing the "contig accession", "start position", "end position"
5. Filtered *T. baileyi* contigs (sense strand) to only include those with both restriction enzyme recognition sites.
6. Filtered *T. baileyi* contigs (antisense strand) to only include those with both restriction enzyme recognition sites.
7. For the set of contigs containing both recognition sites, extract the region between each pairwise combination of RE sites.
8. Filter extracted regions to keep only those with length between 900–1000bp.
9. BLAST (tblastx, tblastn, blastx, blastn?) each sequence in the set of 900-1000bp extracted regions to search within each snake genome
10. Keep the set of single-copy sequences present in all snakes genomes, and design probes for these target loci.

Most of the important files and scripts for selecting ddRAD-like loci are in the zip file **RandomLoci.zip**

**GREPmethod_RandomLociSelection.txt**

**GREPmethod_RandomLociSelection_Part1_MakeHitTable_cluster.R** (run using **GREPmethod_RandomLociSelection_Part1_MakeHitTable_cluster.sh**)

**GREPmethod_RandomLociSelection_Part2_MakeProposedLociTables_cluster.R** (run using **GREPmethod_RandomLociSelection_Part2_MakeProposedLociTables_cluster.sh**)

Proposed target loci using PstI and HpaII recognition sites (didn't use any of these): **proposedLoci_CTGCAG-CCGG_output.txt** (result: 660,088 possible targets).

Proposed target loci using SbfI and EcoR1 recognition sites (USED 328 of these): **proposedLoci_CCTGCAGG-GAATTC_output.txt** (result: 1,178 possible targets); These were reduced to 1,170 possible targets: **compProposedLoci_CCTGCAGG-GAATTC_output.txt**.

Set of 900–1000bp regions of the Sense Strand containing Sbfi and EcoRI recognition sites: **ddRAD-like-loci_SenseStrand_SbfI-EcoRI_900to1000bp_PASSED_HitTable.txt**

<a name="Methods.SelectingMHC"></a>
## Selecting MHC loci:

I used grep to search within the annotation table of the *T. sirtalis* genome (**ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3**) for CDS features of major histocompatibility genes, using the following grep search terms: (1) "MHC", (2) "major histocompatibility". Results = 86 CDS regions corresponding to exons of 19 genes (**ref_Thamnophis_sirtalis-6.0_top_level_JLW_immune-loci-CDS.gff3**).

Target loci included the entire CDS region plus approximately equal amounts of upstream and downstream non-coding DNA (0-60nt). The amount of flanking non-coding DNA targetted was determined by the size of baits (120nt baits). Specifically, the following R script was used to (1) define coordinates of target loci, and (2) to extract MHC targets from the *T. sirtalis* genome and save them to the file **MHC-target-loci_preliminary.fasta**:

```
library(RCurl)
source("~/SnakeCap_functions.R")

annotations.MHC  <- read.table(file="~/ref_Thamnophis_sirtalis-6.0_top_level_JLW_immune-loci-CDS.gff3",sep="\t")
MHC.contigs      <- annotations.MHC[,1]
CDS.start        <- apply(X=annotations.MHC[,4:5],MARGIN=1,FUN=min) ### for each row, the start is position is defined as the smaller of the two values in column 4 and 5
CDS.end          <- apply(X=annotations.MHC[,4:5],MARGIN=1,FUN=max) ### for each row, the end position is defined as the larger of the two values in column 4 and 5
CDS.length       <- (CDS.end-(CDS.start-1))                         ### lengths of CDS regions

bait.length      <- 120 ### nucleotide length of baits

target.length    <- ceiling((ceiling(CDS.length/(bait.length))*(bait.length)))+1
target.start     <- ceiling(CDS.start-(target.length-CDS.length)/2)
target.end       <- ceiling(CDS.end+(target.length-CDS.length)/2)

MHC.targets <- get_ncbi_sequences(outfile="~/MHC-target-loci_preliminary.fasta",
        accessionList = MHC.contigs,
	startList = target.start,
	endList = target.end,
	rettype = "fasta",
)
```

Next, I used blastn to search for the MHC target loci in the *T. sirtalis* genome (results in hit table **MHC_86-targets-vs-Thamnophis_HitTable.csv**).

Examining the hit tables in R:

```
# read in the blastn hit table
mhc.hits.Thamnophis            <- read.csv(file="~/MHC_86-targets-vs-Thamnophis_HitTable.csv",header=F)
colnames(mhc.hits.Thamnophis)  <- c("query.acc.ver", "subject.acc.ver", "percent.identity", "alignment.length", "mismatches", "gap.opens", "query.start", "query.end", "subject.start", "subject.end", "evalue", "bit score")

# read target MHC sequences
library(Biostrings)
mhc.dna        <- readDNAStringSet(filepath="/Users/alyssaleinweber/Documents/Jeff_SequenceCapture-GitHub/MHC-target-loci_preliminary.fasta")
names(mhc.dna) <- gsub(pattern=" .*",replacement="",names(mhc.dna))

# create a vector for query length (mhc target length) corresponding to the rows of the hit table so that query coverage can be calculated
queryLength <- vector(mode="numeric",length=nrow(mhc.hits.Thamnophis))
for(i in 1:length(mhc.dna)){
	queryLength[which(mhc.hits.Thamnophis[,1]==names(mhc.dna[i]))] <- width(mhc.dna[i])
}
# calculate query coverage for each blastn match
queryCoverage <- round((mhc.hits.Thamnophis[,"query.end"]-(mhc.hits.Thamnophis[,"query.start"]-1))/queryLength,digits=4)

## Summarize the blastn hit table
# filter hit table to only include matches with percent.identity and queryCoverage > 90%
# mhc.hits.filtered.90 <- mhc.hits.Thamnophis[which(mhc.hits.Thamnophis[,"percent.identity"] > 90 & queryCoverage > 0.9),]

x.vals      <- c(0,60,70,80,90,95,100)     ### x = threshold for percent.identity and percent queryCoverage
n.matches   <- c(0:100)                    ### range for number of hits in genome to check
summary.mat <- matrix(data="",ncol=length(x.vals),nrow=length(n.matches))

for(i in 1:length(x.vals)){
	mhc.hits.temp <- mhc.hits.Thamnophis[which(mhc.hits.Thamnophis[,"percent.identity"] >= x.vals[i] & queryCoverage >= (x.vals[i]/100)),]
	for(j in 1:length(n.matches)){
		summary.mat[j,i] <- length(which(table(mhc.hits.temp[,1])==n.matches[j]))
	}
}

length(which(table(mhc.hits.filtered.100[,1])==0))


```

Summary of the MHC filtered hit table. The table below shows the number MHC loci (columns 2-5) with at least n matches (column 1) in the genome of *T. sirtalis* with at least *x* percent coverage and percent identical sites across the covered region.

  n matches | MHC loci: x=0 | x=60 | x=70 | x=80 | x=90 | x=95 | x=100
 ---|---|---|---|---|---|---|---
  1    | 28   | 30   | 31   | 35   | 50   | 58   | 70
  2    | 27   | 26   | 27   | 30   | 31   | 27   | 15
  3    |  4   |  4   |  6   |  5   |  2   |  0   |  0
  4    |  2   |  5   |  3   |  2   |  2   |  0   |  0
  5    |  6   |  6   |  6   |  2   |  0   |  0   |  0
 &gt; 5| 19   | 15   | 13   | 12   |  0   |  0   |  0

The MHC loci with genomic coordinates NW_013661433.1:30811-30931 and NW_013659533.1:83942-84062 were dropped from the set of potential target loci because contained very short coding regions (4bp and 19bp, respectively).

The remaining 84 MHC target loci were submitted to Arbor Biosciences for ultrastringent filtration and bait design (WeinellEntry1815-1898):
 - 29 MHC loci failed ultrastringent filtration, and therefore no baits designed for these loci (**Version1_ZeroBaitCoverageLoci.tsv**)
 - 16 other MHC loci were filtered because their baits were all non-specific within the genomes of *T. sirtalis* (**blast results files: XXXXXX)**
 - 12 other MHC loci were already picked as targets (as REEs), and therefore I removed these duplicated targets (**Version1_removed-loci_duplicate-targets.tsv**). 
 
The remaining 27 MHC loci (non-REEs) and six others that were targetted as REEs (entries 248, 559, 728, 787, 891, and 1944) were included in the final set of target loci for which baits were synthesized.

<!--
Version1_Target-loci_Jeff-Weinell_10Sep2018.fasta: 1802 REEs (non-MHC), 12 REE/MHCs, 84 MHCs
Version1_ZeroBaitCoverageLoci.tsv: 70 REEs (non-MHC), 29 MHC loci
Version1_removed-loci_baits-nonspecific.tsv: 8 REEs (non-MHC), 2 MHC loci
Version1_removed-loci_duplicate-targets.tsv: 12 MHC loci
Version2_additional-targets_20Sep2018.txt: 337 REEs (non-MHC), 1 REE (MHC)
Version2-loci-removed_baits-nonspecific.tsv: 111 REEs (non-MHC), 4 REE-MHCs, 14 MHC loci, and 12 UCEs
Version3_additional-targets_Entry3153to5735_4Oct2018.fasta: 2,025 short exons, 98 scalation loci, 132 vision loci, 328 ddRAD-like loci
Version3-loci-removed_ZeroBaitCoverageLoci.tsv: 67 short exon fragments, 7 vision loci
Version3-loci-removed_others.tsv: 209 REEs (non-MHC), 3 MHC-REEs, 1,958 short exon fragments, 10 UCEs, 3 scalation, and 6 vision loci.
See the README file in ArborFiles folder for a description about how the bait kit changed after working with Arbor.
-->

<a name="Methods.SelectingScalation"></a>
## Selecting scalation loci:

I targeted a subset of the genes included in the study by Holthaus et al. (2017). In that study, the authors identified homologous genes of the Epidermal Differentiation Complex (which are putatively involved in scalation) of *Python bivittatus* and *Ophiophagus hannah*. I downloaded the *Ophiophagus* scalation gene sequences using the table of genomic coordinates provided by Holthaus et al. (2017), and then used tblastn to search for and obtain homologous loci in *T. sirtalis*, *Protobothrops mucrosquamatus*, and *Crotalus horridus*.


<a name="Methods.SelectingVision"></a>
## Selecting vision loci:

I used blastn to search for the vision loci probes from Schott et al. (2017) (which were from *Anolis*, *Columba*, *Gallus*, and *Pelodiscus*, *Sceloporus*, or *Python*) within the snake genomes. Most of the SnakeCap probes for these loci are designed from *Ophiophagus* (n = 88), but some probes were designed from *Thamnophis* (n = 21), *Protobothrops* (n = 5), *Pantherophis* (n = 3), or *Python* (n = 2), when blastn of Schott et al 2017 probes did not yield a strong match in *Ophiophagus*.

<a name="ultrastringentFiltering"></a>
<a name="ProbeSynthesis"></a>
## Ultra-stringent filtration and probe synthesis

After choosing the target loci, probes were designed by Arbor Biosciences with the following specifications: 50% tiling, 120nt/probe; 20,020 probes in total. See **Target-loci_Coverage_graph_22October2020.pdf** for a visual summary of target loci, probes, probe coverage, and features of loci including genes, mRNA/transcribed regions, and protein-coding (CDS) regions. This graph was generated with **graph_target_and_features.R** and then filesize reduction in Adobe Acrobat.

**Table X.** Genomes from which synthesized baits were designed from.

species | genome assembly accession | contig prefixes | n target loci with baits designed from species | n bait-containing contigs
----|----|----|----|----
*Thamnophis sirtalis* | GCF_001077635.1 | NW_0136; LFLD01 | 2,619  | 1,412 
*Ophiophagus hannah* | GCA_000516915.1  | AZIM01 | 151 | 86 
*Crotalus horridus* | GCA_001625485.1 | LVCR01 | 2  | 2 
*Python bivittatus* | GCF_000186305.1 | NW_0065; AEQU02 | 18  | 7 
*Probothrops mucrosquamatus* | GCF_001527695.1 | NW_0153; BCNE02 | 9  | 9
*Pantherophis guttatus* | GCA_001185365.1 | JTLQ01 | 2  | 2
*Thermophis baileyi* | GCA_003457575.1 | QLTV01 | 328  | 229

<a name="Sampling"></a>
## Taxa sampled

The following species were sequenced: *Achalinus spinalis*, *Aparallactus capensis*, *Aplopeltura boa*, *Elapsoidea sundevalli*, *Boiga irregularis*, *Buhoma depressiceps*, *Cerberus schneideri*, *Chamaelycus fasciatus*, *Cyclocorus lineatus*, *Oxyrhabdium* cf. *modestum*, *Oxyrhabdium leporinum*, *Oxyrhabdium modestum*, *Prosymna visseri*, *Psammodynastes pulverulentus*, *Pseudaspis cana*, *Pseudoxyrhopus tritaeniatus*, *Rhamphiophis oxyrhynchus*, *Scolecophis atrocinctus*, *Tantilla taeniata*.

<a name="LibraryPrep"></a>
## Sequence capture library prep

Conducted by Arbor Biosciences; eight samples/pool; ...

<a name="DNASequencing"></a>
## DNA sequencing

Novogene Illumina HiSeqX; paired-end sequencing, read length 150nt, insert size 400nt?.

<a name="PostSequencing"></a>
## Post-sequencing

<a name="Demultiplexing"></a>
#### Demultiplexing

To demultiplex raw DNA sequence reads from Illumina, I either used used **bcl2fastq** (Illumina; not sure if I actually used this)... or (more likely) the **bbmap** script **demuxbyname.sh** (Bushnell, 2014).

Running **demuxbyname.sh**:

```
demuxbyname.sh in="HF10_N2_USPD16097067_HY25JBBXX_L4_1.fq" in2="HF10_N2_USPD16097067_HY25JBBXX_L4_2.fq.gz" out=out_%_1.fq out2=out_%_2.fq prefixmode=f names=GGTTACTG+CGAACTTA,TGCGTACA+CATTCGCT,TAGCGTTG+GGTGTTCG,ACTTCCTG+GGTGTCCG,CCGAAGTA+ACCGGTTC,TCGAAGCT+TGAACAGG,GAATCTGG+TTCTGGTG,AATCGGCG+TGAAGCCA,GGCTACAG+ATCACTAC,AAACATCG+TATGCGGT,GGACTATT+ATTCAGAA,TCTACGAC+TGCCTATG,GTGTTGTA+CGTGGACA,TAATGCGC+CAGTGTGG,CGTCAACG+GAGATTCC,TTTCATAG+GCCAGAGG,CTTGCTAT+AAGTCTCC,TGGAGCTG+AAGGATAA,TGCCTATG+GCCAAGAC,CACCACGG+TTGTTCTC,TTGTAGAT+TAAGATTA,GTTCATCT+CAATGATG,AACAGTTG+ACTCCTCC,ATTCAGAA+GGAATTAA,GTGCTTAT+GATCTGCC,GAATCAAT+AGCATCAG,ACCGTAGT+TCGTAGAT,ATCCGCAG+ACCTCCAA,CGTACGTT+TCTCAATT,TTGCCATC+AGCATATT,AGGTGGTC+ACGCTTAT,ACCTTAGA+CTTTTTGA,CATGATGA+TACCACCA,ACAGCAGA+ACTACTTA,TGAAGAGA+TTCGTTCT,ACGGCCGC+TAATGCGC,GATGCCGG+TACAGGTC,AACACATA+CTACGCAT,GTTATATA+GGTTACTG,GCGCCGTG+ATACTACT,TTCGAACC+TCGCCGGC,ACCTCCAA+TACTGTTA,ACTCCTCC+AGTCTTCT,TAACATAG+TGAGTTAG,TATGCGGT+ACTAGCTC,AAGACTGT+GAATCTGT,AAGGTAGG+GGCGAGGA,CCGTGAGA+CAATCGAA,ACGCTTAT+CTACCTTG,ACCGCTAT+TAGCACTT,TGTCTAAC+TGAATGCG,GAATATCC+ATTCAGCG,GAATTCGT+TAAGTACC,GGATTAGG+CGGCGTAA,TTAGAGTC+TGTTAGAC,GAAGTCTT+CAGAAGAT,AAGACGAA+ACACACCT,GACTAGTA+AAGTACAG,GTTCTACT+TGACTACT,TTGCGTAC+TAGCCGAT,CTGTCGAG+GATAGACA,GATAGAGG+AGTTACAT,CTGGAAGC+GAGCTGAA,TAGGATGA+GCTGCATG,GGCGAGGA+TTGGCAGG,ATAGCGAC+CGTCCGTG,AAGTCTCC+CAACTGCT,AATGTTGC+AAGGCAAT,GGTGTCCG+GTGGCTAC,AGCAGGAA+CTCTTGAA,TTATCTAC+GGATAATA,AATCGTTA+TGTTCTCC,GGTTGACG+TGCGTGAA,TGTGGTTG+GTCTGTGC,TAGCCGAT+TTAGGTTG,GTGATTCC+GCATAATT,ATCACAGA+AAGTTATC,CCTAGCCA+AGGTGGTC,ACCATTAA+AACGTGTA,CTACGAGG+CGGATTGC,GAGCTGAA+TAACATAG,AATGCAAA+ATCCGCAG,CTGGCCTC+AAGTCGTG,TAGCACTT+GGCCATCA,GCCAGAGG+GATCTCTT,CGTATTGG+TGCTTGTC,CAGAAGAT+GGTTGACG,GCAGCCTC+CCGCTACA,TGGCTCAG+GCCTGTTC,TGGTCATT+GGACTCTG,GTCGCTGT+GCAGCATA,CAAGCCGC+ACCATTAA,CGCGCCAA+GTAAGGTG,GAAGAGGC+ACAGCCTT,ACCTGACT+ACCGCTAT,CTACTGAC+ACCATATC,CGGATAAC+GAGTTAGC,ACTCTACG+GTCCACTC,TCGGATGT+TGAAGAAT,TCTATCAG+CAGGTTCC,TTGTGTTC+AGGCTATA,AAGGATAA+GCATGGCT,CGCATACA+ACCTTAGA,TTCAAGAA+TTCGGCCG,ACAGGCAG+TAGCTTGT,TTGGTGGC+ATTACTCG,CGGTGGTA+TTATGTAT,GCACTTGG+ATCATTCC,ACGTAGTC+ATTGGCTC,AAGAGAGC+TGAGGCGC,CAGATATT+TATAAGTC,AATCCGTT+AACACATA,AGTGTGTC+AAACATCG
```
<a name="ProcessingReads"></a>
#### Processing sequence reads

To processes sequence reads (assemble contigs for each sample) I followed the FrogCap pipeline (Hutter et al., 2019), which involved running the following R scripts:
- [01_Pre_Process_Reads_Apr10.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/01_Pre_Process_Reads_Apr10.R?token=AJJOG2RAQLT3L3SZXBRDB5273N22S)
- [02_Assemble_Spades_Apr18.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/02_Assemble_Spades_Apr18.R?token=AJJOG2TYC7VHJWH7EG2MADC73N26O)
- [03_Target-loci_matching_20Feb2020.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/03_Target-loci_matching_20Feb2020.R?token=AJJOG2QATDXEF67C6ZX3CPK73N27G) (= **03_Probe-Matching.R** of Hutter et al., 2019) Summary of results:[Sample-Assessment.tsv](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/Sample-Assessment.tsv?token=AJJOG2XJAXWKGTBXVUBWCCK73N2LQ)
- [03-2_Data-subsetting_JLW.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/03-2_Data-subsetting_JLW.R?token=AJJOG2RI6Q7NRMHJCCYPRM273N3AM) (This is an extra step not in Hutter et al., 2019)
- [04_Loci_alignment_1May2019.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/04_Loci_alignment_Aug3.r?token=AJJOG2VU4DTHM26XFKK5QHC73N3FW)
- [05_mtgenome_assembly_May8.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/05_mtgenome_assembly_May8.R?token=AJJOG2ULM27MZDC5CDRD7DS73N3GM)
- [06_Trim_Align_Aug14.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/06_Trim_Align_Aug14.R?token=AJJOG2T33LH4U5ZOMF27OZ273N3G6)
- [07_Concat_CompleteMatrix_Aug29.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/07_Concat_CompleteMatrix_Aug29.R?token=AJJOG2WK35XFGXC2XA2JEWS73N3IG)
- [07-2_IQTREE_1May2019.R](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/PostSequencing/07-2_IQTREE_1May2019.R?token=AJJOG2XRMBXFMEFDSBULZG273N3I2)
- **Sort_Alignments_by_LocusType.R**

<a name="DNA.Alignment"></a>
#### DNA alignment

**Generating an unpartitioned, multiple-species alignment for each captured locus:**

I used the MAFFT algorithm for multiple alignment, as implemented in the R package rMSA

```
mafft(final.locus,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
```

**Generating partitioned and subset alignments:**

I sorted each of the unpartitioned target+flanking alignments into one of six directories by locus type (REEs, UCEs, ddRAD-like, Immune, Scalation, Vision), and then ran an R script to generate nine types of alignments per locus, which differ in how data are subset or partitioned; the nine alignments are:

- CDS_only: only includes protein-coding regions
- CDS_FirstCodonPosition: only
- CDS_SecondCodonPosition
- CDS_ThirdCodonPosition
- AminoAcids (an alignment of amino acid sequences from translating the CDS_only alignment)
- Upstream_noncoding   (only includes noncoding DNA upstream of the first captured exon)
- Downstream_noncoding (only includes noncoding DNA downstream of the last captured exon)
- All_noncoding (includes only the non-CDS regions)
- All_parts (includes all data types, and a partition file is created)

*Note*: Most of the UCEs and ddRAD-like loci are non-protein-coding, and therefore CDS only subset alignments were not created for these loci; 87/907 UCEs and 90/328 ddRAD-like loci contain CDS regions. In contrast, 91/119 vision loci, 92/95 scalation loci, 27/27 immune loci and 1653/1653 REEs contain CDS regions.

In R, load in SnakeCap functions and required packages:

```
source("~/SnakeCap_Functions.R") # load in the functions in this file
library(Biostrings)
library(stringr)
library(ape)
library(data.table)
```
Run make.partitioned.alignment separately for REEs, UCEs, ddRAD-like, MHC, scalation, and vision genes:

```
### REEs
make.partitioned.alignment(InputAlignmentFolder="~/WholeExon/unpartitioned/", output.dir="~/WholeExon/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")

### UCEs
make.partitioned.alignment(InputAlignmentFolder="~/UCE/unpartitioned/", output.dir="~/UCE/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")

### ddRAD-like loci
make.partitioned.alignment(InputAlignmentFolder="~/ddRAD-like/unpartitioned/", output.dir="~/ddRAD-like/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")

### Immune (MHC) genes:
make.partitioned.alignment(InputAlignmentFolder="~/Immune/unpartitioned/", output.dir="~/Immune/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")

### Scalation genes:
make.partitioned.alignment(InputAlignmentFolder="~/Scalation/unpartitioned/", output.dir="~/Scalation/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")

### Vision genes:
make.partitioned.alignment(InputAlignmentFolder="~/Vision/unpartitioned/", output.dir="~/Vision/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")
```

<a name="PhylogeneticAnalyses"></a>
#### Phylogenetic and Network Analyses

- IQTREE...
- ASTRAL III...
- Icytree...
- PhyloNet...
- NetView...
- Dendroscope...

<a name="Results"></a>
## Results

<a name="References"></a>
## References:

SPADES...

Bushnell, B. 2014. BBMap: A Fast, Accurate, Splice-Aware Aligner...

Holthaus K.B., Mlitz V., Strasser B., Tschachler E., Alibardi L., and L. Eckhart. 2017. Identification and comparative analysis of the epidermal differentiation complex in snakes. *Scientific Reports* 7, 45338. doi: http://doi.org/10.1038/srep45338.

Hutter C.R., Cobb K.A., Portik D., Travers S., Wood Jr. P.L., and R.M. Brown. 2019. FrogCap: A modular sequence capture probe set for phylogenomics and population genetics for Anurans, assessed across multiple phylogenetic scales. *bioRxiv* 825307. doi: https://doi.org/10.1101/825307.

Schott R.K., Panesar B., Card D.C., Preston M., Castoe T.A., and B.S.W. Chang. 2017. Targeted Capture of Complete Coding Regions across Divergent Species. *Genome Biology and Evolution* 9(2), 398–414. doi: http://doi.org/10.1093/gbe/evx005

Streicher J.W., and J.J Wiens J.J. 2017. Phylogenomic analyses of more than 4000 nuclear loci resolve the origin of snakes among lizard families. *Biology Letters* 13, 20170393. doi: http://doi.org/10.1098/rsbl.2017.0393.











