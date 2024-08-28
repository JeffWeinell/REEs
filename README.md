# SnakeCap Sequence Capture Probe Set

# Contents

[Description of SnakeCap probe set](#Description)

[Selecting loci to capture](#Methods.SelectingTargetLoci)
  - [Rapidly-evolving Exons (REEs)](#Methods.SelectingREEs)
  - [Ultraconserved-elements (UCEs)](#Methods.SelectingUCEs)
  - [ddRAD-like loci](#Methods.SelectingddRAD)
  - [Major Histocompatibility loci (MHCs)](#Methods.SelectingMHC)
  - [Scalation loci](#Methods.SelectingScalation)
  - [Vision loci](#Methods.SelectingVision)

[Ultrastringent filtration and Probe Synthesis](#ProbeSynthesis)

[References](#References)

<a name="Description"></a>
# Description of SnakeCap probe set

The probe set includes 20,020 probes for 3,128 single-copy loci (1,517,011 nt) in genomes of a diverse group of snakes.

Target loci include rapidly evolving exons (REEs), ultra-conserved elements (UCEs), ddRAD-like loci, major histocompatibility complex (MHC) genes, vision-associated genes, and scalation-associated genes.

REEs include one or more entire exons, one or both exon-flanking regions, and range in length from 121 to 7,501 nt. I used a modified version of the FrogCap pipeline (Hutter et al., 2019) to select the optimal set of REEs from an alignment of snake exomes.

SnakeCap UCEs are a subset of the *Micrurus fulvius* UCEs from Streicher and Wiens (2017).

ddRAD-like loci are shared, single-copy loci identified from in-silico ddRAD using recognition sites for SbfI and EcoRI restriction enzymes.

MHC, vision, and scalation loci included entire or subregions of genes. Putative gene functions are from results of previous studies.

**Table 1**. For each type of locus: genomic region targeted, number of loci (nloci), number of nucleotides targeted (nt), and nucleotide lengths (nt/locus) of the shortest and longest loci.

Locus type | Region targeted | nloci | nnt | nt/locus (min–max)
---- | ---- | ---- | ---- | ----
Rapidly evolving exons (REEs) | Usually the entire exon + 0–60nt each of 5' & 3' flanking regions. | 1,652 | 996,369 | 121–7,501
Ultra-conserved elements (UCEs) | Entire UCE region previously identified | 907 | 143,075 | 120–161
ddRAD-like | in silico ddRAD selected loci | 328 | 271,505 | 120–996
MHC | Exons + 0–60nt each of 5' & 3' flanking regions of major histocompatibility complex (MHC) genes. | 27 | 5,354 | 121–364
Vision | 160nt, including ≤ 70nt upstream of start codon if targeted exon is first exon. Usually, first 160nt of exon targeted. | 119 | 18,857 | 120–170
Scalation | 1100nt, including 1000nt of promoter region + first 100nt of first exon. | 95 | 81,851 | 125–1,101
All loci |  | 3,128 | 1,517,011 | 120–7,501 (mean = 531.62)  


<a name="Methods.SelectingTargetLoci"></a>
# Selecting loci to capture

<a name="Methods.SelectingREEs"></a>
## Selecting the set of target REEs

#### Overview: 

I used the following R functions (shown as packageName::functionName):

REEs::load.gff --> REEs::filter.gff --> REEs::get.seqs.from.gff --> Biostrings::writeXStringSet --> REEs::blast --> REEs::reportBestMatches --> REEs::get.seqs.from.blastTable --> REEs::makeStatsTable --> REEs::pick.loci --> (then functions in step 8  to get REEs + small region of exon-flanking DNA).

#### Details:

<!--
I downloaded the [*Thamnophis sirtalis* genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz) and its associated annotation table: [GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz) (n = 559,130 features annotated). Then, I renamed the contigs in the genome file to have the following format: **Thamnophis_sirtalis_GCF_001077635.1_read1**, **Thamnophis_sirtalis_GCF_001077635.1_read2**, etc., and saved this renamed genome in sequential fasta format: [ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3.zip](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/exomes/ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3.zip). The two-column, tab-delimited table [Scaffold-Name-Key.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/exomes/Scaffold-Name-Key.txt?token=AJJOG2UQ6MDA7UY2U4R6BFS7ZDZYS) includes the new contig name in the first column and the original contig name in the second column:
-->

Step-by-step description of pipeline used to select target REEs:

**1**. I used the R functions REEs::load.gff and REEs::filter.gff to filter the *Thamnophis sirtalis* genome feature table [GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz) to include only CDS features ≥ 120bp.

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

**2**. I used the function REEs::get.seqs.from.gff to extract the sequences corresponding the CDS features in the table generated in step 1. 

```
### URL to the Thamnophis sirtalis genome (fasta formatted sequences).
Thamnophis.sirtalis_genome.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz"

### Extract the sequences for the loci in the filtered GFF (Thamnophis.sirtalis_GFF_CDS_longer120bp from step 1).
Thamnophis.sirtalis_exome   <- REEs::get.seqs.from.gff(input.seqs=Thamnophis.sirtalis_genome.path,input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp) 

### Save extracted sequences
writeXStringSet(x=Thamnophis.sirtalis_exome,filepath="./Thamnophis_sirtalis_exome_longer120bp.fas")

```
Output DNA sequences were saved in fasta format and can be downloaded here: [Thamnophis_sirtalis_exome_longer120bp.fas](https://osf.io/v7ecb/download).

**3**. I used TBLASTX as implemented in the REEs::blast function to search for each *T. sirtalis* CDS sequence (from step 2) in each squamate genome listed in Table 2. I saved up to 50 matches per query sequence.

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
Vipera.berus.genome_url                 <- REEs::datasets(1)[which(datasets(1)[,1]=="Vipera berus"),2]
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
Vipera.berus.50hits                 <- REEs::blast(method="tblastx",subject=Vipera.berus.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Vipera.berus.tblastx.exons.50hits.txt")
Thamnophis.sirtalis.50hits          <- REEs::blast(method="tblastx",subject=Thamnophis.sirtalis.genome_url,query=Thamnophis.sirtalis_exome,table.out="./Thamnophis.sirtalis.tblastx.exons.50hits.txt")
```
The TBLASTX output tables can be downloaded here (these are xz compressed): [Anolis.carolinensis.tblastx.exons.50hits.txt.xz](https://osf.io/3a4be/download), [Gekko.japonicus.tblastx.exons.50hits.txt.xz](https://osf.io/u2ec5/download), [Pogona.vitticeps.tblastx.exons.50hits.txt.xz](https://osf.io/9ehyf/download), [Crotalus.horridus.tblastx.exons.50hits.txt.xz](https://osf.io/j4gvk/download), [Crotalus.mitchellii.tblastx.exons.50hits.txt.xz](https://osf.io/5vxt7/download), [Ophiophagus.hannah.tblastx.exons.50hits.txt.xz](https://osf.io/jr5zd/download), [Pantherophis.guttatus.tblastx.exons.50hits.txt.xz](https://osf.io/tu3ve/download), [Protobothrops.mucrosquamatus.tblastx.exons.50hits.txt.xz](https://osf.io/7xmr4/download), [Python.bivittatus.tblastx.exons.50hits.txt.xz](https://osf.io/8s5na/download), [Vipera.berus.tblastx.exons.50hits.txt.xz](https://osf.io/48x7e/download), [Thamnophis.sirtalis.tblastx.exons.50hits.txt.xz](https://osf.io/ctz73/download).

**4**. I used the REEs::reportBestMatches function to filter the hit tables generated by TBLASTX to include only the best match (max bitscore) per query. Matches with bitscore < 50 were dropped. If bitscores of the best and second-best matches differed by < **5**, then the query locus was filtered.

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
best.hits.Vipera.berus                  <- REEs::reportBestMatches(input.table=Vipera.berus.50hits,output.table.path="./Vipera.berus.best.hits.txt")
```

Output tables containing the set of best matches: [Anolis.carolinensis.exons.best.hits.txt](https://osf.io/a5b8n/download), [Gekko.japonicus.exons.best.hits.txt](https://osf.io/a9p23/download), [Pogona.vitticeps.exons.best.hits.txt](https://osf.io/xrycw/download), [Crotalus.horridus.exons.best.hits.txt](https://osf.io/zrndt/download), [Crotalus.mitchellii.exons.best.hits.txt](https://osf.io/bhnwy/download), [Ophiophagus.hannah.exons.best.hits.txt](https://osf.io/y49cv/download), [Pantherophis.guttatus.exons.best.hits.txt](https://osf.io/2vw9f/download), [Protobothrops.mucrosquamatus.exons.best.hits.txt](https://osf.io/aqps2/download), [Python.bivittatus.exons.best.hits.txt](https://osf.io/wxthd/download), [Vipera.berus.exons.best.hits.txt](https://osf.io/en4yg/download),[Thamnophis.sirtalis.exons.best.hits.txt](https://osf.io/mekqz/download).

**5**. I used get.seqs.from.blastTable function (REEs package) to extract sequences for the best matches identified in step 4.

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
Vipera.berus.best.hits.seqs                 <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Vipera.berus,input.seqs=Vipera.berus.genome_url,output.path="./Vipera.berus.tblastx.best.hits_seqs.fas")
Thamnophis.sirtalis.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=best.hits.Thamnophis.sirtalis,input.seqs=Thamnophis.sirtalis.genome_url,output.path="./Thamnophis.sirtalis.tblastx.best.hits_seqs.fas")
```

The output sequences were saved in fasta format and can be downloaded here: [Anolis.carolinensis.tblastx.best.hits_seqs.fas](https://osf.io/mabyw/download), [Gekko.japonicus.tblastx.best.hits_seqs.fas](https://osf.io/8etnk/download), [Pogona.vitticeps.tblastx.best.hits_seqs.fas](https://osf.io/cxey4/download), [Crotalus.horridus.tblastx.best.hits_seqs.fas](https://osf.io/apvx2/download), [Crotalus.mitchellii.tblastx.best.hits_seqs.fas](https://osf.io/tm4jq/download), [Ophiophagus.hannah.tblastx.best.hits_seqs.fas](https://osf.io/5d7wa/download), [Pantherophis.guttatus.tblastx.best.hits_seqs.fas](https://osf.io/byj3d/download), [Protobothrops.mucrosquamatus.tblastx.best.hits_seqs.fas](https://osf.io/6yk2b/download), [Python.bivittatus.tblastx.best.hits_seqs.fas](https://osf.io/ywdq9/download), [Vipera.berus.tblastx.best.hits_seqs.fas](https://osf.io/8vpq5/download), [Thamnophis.sirtalis.tblastx.best.hits_seqs.fas](https://osf.io/g8rqd/download).

**6**. I used the makeStatsTable function (REEs package) to perform multiple sequence alignment (MAFFT algorithm) for the best matches to each query exon sequence and to create a table of alignment statistics for each alignment. The information and statistics in the output table are summarized in Table 3.

```
### Create a character vector holding the paths to each of the ".*.tblastx.best.hits_seqs.fas" files generated in step 5. 
input.seqs.paths <- c("Anolis.carolinensis.tblastx.best.hits_seqs.fas", "Gekko.japonicus.tblastx.best.hits_seqs.fas", "Pogona.vitticeps.tblastx.best.hits_seqs.fas", "Crotalus.horridus.tblastx.best.hits_seqs.fas", "Crotalus.mitchellii.tblastx.best.hits_seqs.fas", "Ophiophagus.hannah.tblastx.best.hits_seqs.fas", "Pantherophis.guttatus.tblastx.best.hits_seqs.fas", "Protobothrops.mucrosquamatus.tblastx.best.hits_seqs.fas", "Python.bivittatus.tblastx.best.hits_seqs.fas", "Vipera.berus.tblastx.best.hits_seqs.fas", "Thamnophis.sirtalis.tblastx.best.hits_seqs.fas")

### Align homologous sequences and make a table to summarize data in each alignment (one row per aligned locus).
stats.table.all  <- makeStatsTable(input.seqs=input.seqs.paths, input.gff=Thamnophis.sirtalis_GFF_CDS_longer120bp, output.path="./statsTable_REEs_SnakeCap.txt", species.names=c("Anolis carolinensis","Gekko japonicus","Pogona vitticeps","Crotalus horridus","Crotalus mitchellii","Ophiophagus hannah","Pantherophis guttatus", "Protobothrops mucrosquamatus", "Python bivittatus","Vipera berus","Thamnophis sirtalis"), subgroup=c("Crotalus horridus","Crotalus mitchellii","Ophiophagus hannah","Pantherophis guttatus", "Protobothrops mucrosquamatus", "Python bivittatus","Vipera berus","Thamnophis sirtalis"),alignments.out=alignments.dir, reference.species=11)
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
17|pident.Vipera_berus|Percent of sites identical between *T. sirtalis* query and *V. berus* best match.
18|gene.name|Name of the gene of the *T. sirtalis* query locus. The is the name associated with the query sequence's CDS feature in the *T. sirtalis* GFF genome feature table.
19|locus.length.Thamnophis_sirtalis|Length of the query sequence
20|mean.variable.sites| Mean pairwise fraction of sites variable compared to *T. sirtalis*. This is simply (100-mean.pident)/100
21|min.pident.all|Minimum of percent of sites identical to *T. sirtalis* among the species. In other words, the minimum of columns 8–17. This would be calculated across different columns if different species were used.
22|min.pident.subgroup|Minimum of percent of sites identical to *T. sirtalis* among the snake species compared. In other words, the minimum of columns 11–17. The would be calculated across different columns if the subgroup was defined differently.
23|alignment.width|Width of the alignment. This column was not generated for the SnakeCap dataset because I used an older version of makeStatsTable, which did not include this column in the output table.

**7**. I used the function pick.loci (REEs package) to choose a set of REEs for sequence capture. I applied the following criteria to select REEs: (1) I filtered out loci if minimum percent genetic similarity (mean among snakes) to *T. sirtalis* was < 65% or = 100%; (2) for genes with multiple exons, I kept only the exon with the lowest mean pairwise genetic distance to *T. sirtalis* (16,650 exons pass this step); (3) the maximum number of baits required to target the exons was set to 20,000 (this constraint was imposed by the 20K my-baits kit) with bait size = 120nt and 50% bait tiling; (4) the maximum number of nucleotides to target was set to 1.2Mb (also constrained by the my-baits kit); (5) of the sets of exons that meet criteria 1–4, I chose the set that maximizes the total number of variable sites relative to *T. sirtalis*.

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

**8**. Expand the target region to include noncoding DNA upstream and downstream of the loci identified in step 7. The amount of noncoding DNA included in each expanded target is a function of the bait length (120nt) and exon length, such that the expanded target length is a multiple of the bait length.

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

Six duplicate pairs of REEs (each pair with identical sequences) were present in the output of step 8 (Table 5). Only one of these pairs was recognized/identified (and filtered manually) prior to submitting target sequences to Arbor Biosciences.

The remaining 2,068 REEs (expanded targets) were submitted to Arbor Biosciences for ulstrastringent filtering and probe design, and can be downloaded here: [REEs.expanded.final.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/REEs/REEs.expanded.final.fas). These were actually submitted to Arbor in two batches: [Version1_Target-loci_Jeff-Weinell_10Sep2018.fasta](https://github.com/JeffWeinell/SnakeCap/raw/main/ArborFiles/Version1_Target-Loci_Jeff-Weinell_10Sep2018.fasta) and [Version2_additional-targets_20Sep2018.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/ArborFiles/Version2_additional-targets_Entry1899to3152_20Sep2018.fasta). See [ultra-stringent filtering](#ultrastringentFiltering) section.

Arbor performed ultrastringent filtration on the 2,068 REEs retained from step 8 (after removing two identical sequences; pair 1 Table 5). Ultrastringent filtering resulted in the removal of 203 REEs (1,865 REEs retained). Of the 203 REEs that were filtered, 76 were filtered because no baits could be designed for these loci (70 of these are listed in [Version1-loci-removed_ZeroBaitCoverageLoci.tsv](https://git.io/JLiEu) and six are listed in [Version2-loci-removed_ZeroBaitCoverageLoci.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version2-loci-removed_ZeroBaitCoverageLoci.tsv)); and 127 REEs were filtered because all proposed baits were non-specific within the *T. sirtalis* genome (eight of these are listed in [Version1-loci-removed_baits-nonspecific.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version1-loci-removed_nonspecific-baits.tsv) and 119 are listed in [Version2-loci-removed_baits-nonspecific.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version2-loci-removed_baits-nonspecific.tsv)).

Of the 1,865 REEs that passed ultrastringent filtering, 212 were removed to allow a fraction of the 20K baits to be used to target other types of loci (UCEs, MHC genes, scalation genes, vision genes, and ddRAD-like loci), and these removed REEs are listed in the file [Version3-loci-removed_others.tsv](https://github.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version3-loci-removed_others.tsv). Baits for the remaining 1,653 REEs were synthesized by Arbor (mybaits 20K bait kit: product no. 3001160).

**Table 4**. The nine partially expanded targets from output of step 8 and reason why targets were partially rather than fully expanded.
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

**Table 5**. Pairs of REEs having identical sequences that were included in the ouput of the pick.loci function (step 8). Five of these pairs were subsequently filtered, either immediately before or after application of Arbor's ultrastringent filtering algorithm. One pair (pair 6) was not identified until after synthesis and sequencing. The latest version of the pick.loci function has an option to filter REEs if the bitscores of the top matches are too similar according to a user-defined threshold.
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
NW_013661433.1|31138|31498|6|WeinellEntry1857| not filtered
NW_013661433.1|31138|31498|6|WeinellEntry1944| filtered from target reference fasta prior to target matching of assembled loci

<!--
REEs that failed ultrastringent filtration: 70
REEs removed because all baits were non-specific: 123
REEs removed to make room for other loci: 212
Total REEs removed: 70+123+212 = 405
1,653 REEs were synthesized and included in the bait kit.
10 REEs unaccounted for...
-->

<a name="Methods.SelectingUCEs"></a>
## Selecting the set of target UCEs

#### Overview:

Target UCEs include 907 of the 3,260 UCEs previously identified in *Micrurus fulvius* (Streicher and Wiens, 2017; **Table X**). First, I filtered the full set of *Micrurus* UCEs to only include those present in all NCBI snake genomes (n = 2,968 UCEs). Then, I filtered the set of shared UCEs to include only loci in which the *T. sirtalis* sequence included at least 200nt of non-gap sites in the UCE alignment. I sorted the set of passing *T. sirtalis* UCEs by UCE name and then retained only the first 1,000 loci. From each of these 1,000 loci, I extracted the middle 161nt bases to be used at a template for probe design. Arbor Biosciences was able to synthesize probes for 907 of the 1,000 proposed target UCEs.


#### Details:

Step-by-step description of pipeline used to select target UCEs:

**1**. I downloaded the set of *Micurus fulvius* UCEs (n = 3,260) identified by Streicher and Wiens (2017). These were available as a fasta file called [micrurus_UCEs.fa](https://git.io/JLiu2).

**2**. To find *Micrurus* UCEs in each of the other snake genomes, I queried each *Micrurus* UCE against each snake genome using blastn algorithm (saving ≤ 50 matches per query), which was implemented using the R wrapper function REEs::blast. Note: The newest version of the *Pantherophis guttatus* genome (JTLQ00000000.2 as of January 13, 2021) was not available when I originally did these analyses. For reproducibility, "*Pantherophis guttatus* genome" refers to version JTLQ00000000.1. The contigs of version JTLQ00000000.1 are split into three fasta files available from NCBI. I downloaded these files and merged them into one file using the script shown here below. The other genomes can be downloaded as a single file from NCBI.

```
#### Download the three fasta files containing contigs of the Pantherophis guttatus genome (JTLQ00000000.1). These will be saved to your current directory unless you specify otherwise.
utils::download.file(url="https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/JT/LQ/JTLQ01/JTLQ01.1.fsa_nt.gz", destfile="JTLQ01.1.fsa_nt.gz")
utils::download.file(url="https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/JT/LQ/JTLQ01/JTLQ01.2.fsa_nt.gz", destfile="JTLQ01.2.fsa_nt.gz")
utils::download.file(url="https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/JT/LQ/JTLQ01/JTLQ01.3.fsa_nt.gz", destfile="JTLQ01.3.fsa_nt.gz")

### Load the downloaded fasta files into R as DNAStringSet objects.
Pantherophis.guttatus.genome_JTLQ01.1 <- Biostrings::readDNAStringSet("JTLQ01.1.fsa_nt.gz")
Pantherophis.guttatus.genome_JTLQ01.2 <- Biostrings::readDNAStringSet("JTLQ01.2.fsa_nt.gz")
Pantherophis.guttatus.genome_JTLQ01.3 <- Biostrings::readDNAStringSet("JTLQ01.3.fsa_nt.gz")

### Merge the three DNAStringSet objects into one object containing all contigs of the genome.
Pantherophis.guttatus.genome_JTLQ01 <- c(Pantherophis.guttatus.genome_JTLQ01.1,Pantherophis.guttatus.genome_JTLQ01.2,Pantherophis.guttatus.genome_JTLQ01.3)

### Save the genome.
Biostrings::writeXStringSet(x=Pantherophis.guttatus.genome_JTLQ01,filepath="Pantherophis.guttatus_JTLQ01.fsa_nt.gz",compress=TRUE)
```
Result: The *Pantherophis guttatus* genome (merged contigs file) of JTLQ00000000.1 can be downloaded here: [Pantherophis.guttatus_JTLQ01.fsa_nt.gz](https://osf.io/r2xa3/download).

```
### Define URLs to fasta genome sequences
# URL to Pantherophis guttatus genome version JTLQ00000000.1
Pantherophis.guttatus.genome_url <- "https://osf.io/r2xa3/download"
# The URLs to the other genomes are held in a matrix that can be accessed with the REEs::datasets function.
Anolis.carolinensis.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Anolis carolinensis"),2]
Gekko.japonicus.genome_url              <- REEs::datasets(1)[which(datasets(1)[,1]=="Gekko japonicus"),2]
Pogona.vitticeps.genome_url             <- REEs::datasets(1)[which(datasets(1)[,1]=="Pogona vitticeps"),2]
Crotalus.horridus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus horridus"),2]
Crotalus.mitchellii.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus mitchellii"),2]
Ophiophagus.hannah.genome_url           <- REEs::datasets(1)[which(datasets(1)[,1]=="Ophiophagus hannah"),2]
Protobothrops.mucrosquamatus.genome_url <- REEs::datasets(1)[which(datasets(1)[,1]=="Protobothrops mucrosquamatus"),2]
Python.bivittatus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Python bivittatus"),2]
Vipera.berus.genome_url                 <- REEs::datasets(1)[which(datasets(1)[,1]=="Vipera berus"),2]
Thamnophis.sirtalis.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Thamnophis sirtalis"),2]

### Runs BLASTN
Crotalus.horridus.UCEs.50hits            <- REEs::blast(method="blastn",subject=Crotalus.horridus.genome_url, query="./micrurus_UCEs.fa", table.out="./Crotalus.horridus.blastn.UCEs.50hits.txt")
Crotalus.mitchellii.UCEs.50hits          <- REEs::blast(method="blastn",subject=Crotalus.mitchellii.genome_url, query="./micrurus_UCEs.fa", table.out="./Crotalus.mitchellii.blastn.UCEs.50hits.txt")
Ophiophagus.hannah.UCEs.50hits           <- REEs::blast(method="blastn",subject=Ophiophagus.hannah.genome_url, query="./micrurus_UCEs.fa", table.out="./Ophiophagus.hannah.blastn.UCEs.50hits.txt")
Pantherophis.guttatus.UCEs.50hits        <- REEs::blast(method="blastn",subject=Pantherophis.guttatus.genome_url, query="./micrurus_UCEs.fa", table.out="./Pantherophis.guttatus.blastn.UCEs.50hits.txt")
Protobothrops.mucrosquamatus.UCEs.50hits <- REEs::blast(method="blastn",subject=Protobothrops.mucrosquamatus.genome_url, query="./micrurus_UCEs.fa", table.out="./Protobothrops.mucrosquamatus.blastn.UCEs.50hits.txt")
Python.bivittatus.UCEs.50hits            <- REEs::blast(method="blastn",subject=Python.bivittatus.genome_url, query="./micrurus_UCEs.fa", table.out="./Python.bivittatus.blastn.UCEs.50hits.txt")
Vipera.berus.UCEs.50hits                 <- REEs::blast(method="blastn",subject=Vipera.berus.genome_url, query="./micrurus_UCEs.fa", table.out="./Vipera.berus.blastn.UCEs.50hits.txt")
Thamnophis.sirtalis.UCEs.50hits          <- REEs::blast(method="blastn",subject=Thamnophis.sirtalis.genome_url, query="./micrurus_UCEs.fa", table.out="./Thamnophis.sirtalis.blastn.UCEs.50hits.txt")
```

Output tables from BLASTN: [Crotalus.horridus.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.horridus.blastn.UCEs.50hits.txt), [Crotalus.mitchellii.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.mitchellii.blastn.UCEs.50hits.txt), [Ophiophagus.hannah.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus.hannah.blastn.UCEs.50hits.txt), [Pantherophis.guttatus.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis.guttatus.blastn.UCEs.50hits.txt), [Protobothrops.mucrosquamatus.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops.mucrosquamatus.blastn.UCEs.50hits.txt), [Python.bivittatus.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Python.bivittatus.blastn.UCEs.50hits.txt), [Thamnophis.sirtalis.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis.sirtalis.blastn.UCEs.50hits.txt), and [Vipera.berus.blastn.UCEs.50hits.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Vipera.berus.blastn.UCEs.50hits.txt).

**3**. I used the function reportBestMatches to filter the hit tables to include only the best match for each query UCE.

```
UCEs.best.hits.Crotalus.horridus             <- REEs::reportBestMatches(input.table=Crotalus.horridus.UCEs.50hits, output.table.path="Crotalus.horridus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Crotalus.mitchellii           <- REEs::reportBestMatches(input.table=Crotalus.mitchellii.UCEs.50hits, output.table.path="Crotalus.mitchellii.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Ophiophagus.hannah            <- REEs::reportBestMatches(input.table=Ophiophagus.hannah.UCEs.50hits, output.table.path="Ophiophagus.hannah.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Pantherophis.guttatus         <- REEs::reportBestMatches(input.table=Pantherophis.guttatus.UCEs.50hits, output.table.path="Pantherophis.guttatus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Protobothrops.mucrosquamatus  <- REEs::reportBestMatches(input.table=Protobothrops.mucrosquamatus.UCEs.50hits, output.table.path="Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Python.bivittatus             <- REEs::reportBestMatches(input.table=Python.bivittatus.UCEs.50hits, output.table.path="Python.bivittatus.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Thamnophis.sirtalis           <- REEs::reportBestMatches(input.table=Thamnophis.sirtalis.UCEs.50hits, output.table.path="Thamnophis.sirtalis.blastn.UCEs.best.hits.txt")
UCEs.best.hits.Vipera.berus                  <- REEs::reportBestMatches(input.table=Vipera.berus.UCEs.50hits, output.table.path="Vipera.berus.blastn.UCEs.best.hits.txt")

### Getting a list of the shared loci.
shared.UCEs <- intersect.all(list(UCEs.best.hits.Crotalus.horridus$qseqid, UCEs.best.hits.Crotalus.mitchellii$qseqid,UCEs.best.hits.Ophiophagus.hannah$qseqid, UCEs.best.hits.Pantherophis.guttatus$qseqid,UCEs.best.hits.Protobothrops.mucrosquamatus$qseqid, UCEs.best.hits.Python.bivittatus$qseqid, UCEs.best.hits.Vipera.berus$qseqid, UCEs.best.hits.Thamnophis.sirtalis$qseqid))

### Result: 2,168 UCE loci were found in all of the species. The loci missing from one or more species are not filtered here, but are filtered later when the function align.shared.loci is used.

```
Output tables were saved to the files: [Crotalus.horridus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.horridus.blastn.UCEs.best.hits.txt), [Crotalus.mitchellii.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.mitchellii.blastn.UCEs.best.hits.txt), [Ophiophagus.hannah.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus.hannah.blastn.UCEs.best.hits.txt) , [Pantherophis.guttatus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis.guttatus.blastn.UCEs.best.hits.txt), [Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt), [Python.bivittatus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Python.bivittatus.blastn.UCEs.best.hits.txt), [Thamnophis.sirtalis.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis.sirtalis.blastn.UCEs.best.hits.txt), [Vipera.berus.blastn.UCEs.best.hits.txt](https://github.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Vipera.berus.blastn.UCEs.best.hits.txt).


**4**. I used the function REEs::get.seqs.from.blastTable to extract and save the set of best-match UCEs from each genome.
```
#### Read in the tables of best matches generated from step 3 (if these are not already loaded).
UCEs.best.hits.Crotalus.horridus             <- as.data.frame(data.table::fread("./Crotalus.horridus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Crotalus.mitchellii           <- as.data.frame(data.table::fread("./Crotalus.mitchellii.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Ophiophagus.hannah            <- as.data.frame(data.table::fread("./Ophiophagus.hannah.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Pantherophis.guttatus         <- as.data.frame(data.table::fread("./Pantherophis.guttatus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Protobothrops.mucrosquamatus  <- as.data.frame(data.table::fread("./Protobothrops.mucrosquamatus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Python.bivittatus             <- as.data.frame(data.table::fread("./Python.bivittatus.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Thamnophis.sirtalis           <- as.data.frame(data.table::fread("./Thamnophis.sirtalis.blastn.UCEs.best.hits.txt"),header=T)
UCEs.best.hits.Vipera.berus                  <- as.data.frame(data.table::fread("./Vipera.berus.blastn.UCEs.best.hits.txt"),header=T)

#### Define the genome URL paths. The URLs for genomes used in the SnakeCap study can be called using the datasets function of REEs package.
Crotalus.horridus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus horridus"),2]
Crotalus.mitchellii.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Crotalus mitchellii"),2]
Ophiophagus.hannah.genome_url           <- REEs::datasets(1)[which(datasets(1)[,1]=="Ophiophagus hannah"),2]
Pantherophis.guttatus.genome_url        <- REEs::datasets(1)[which(datasets(1)[,1]=="Pantherophis guttatus"),2]
Protobothrops.mucrosquamatus.genome_url <- REEs::datasets(1)[which(datasets(1)[,1]=="Protobothrops mucrosquamatus"),2]
Python.bivittatus.genome_url            <- REEs::datasets(1)[which(datasets(1)[,1]=="Python bivittatus"),2]
Vipera.berus.genome_url                 <- REEs::datasets(1)[which(datasets(1)[,1]=="Vipera berus"),2]
Thamnophis.sirtalis.genome_url          <- REEs::datasets(1)[which(datasets(1)[,1]=="Thamnophis sirtalis"),2]

#### Extracts the sequence of the best match of each UCE from each snake genome.
Crotalus.horridus.best.hits.seqs            <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Crotalus.horridus, input.seqs=Crotalus.horridus.genome_url, output.path="./Crotalus.horridus.UCEs.fas")
Crotalus.mitchellii.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Crotalus.mitchellii, input.seqs=Crotalus.mitchellii.genome_url, output.path="./Crotalus.mitchellii.UCEs.fas")
Ophiophagus.hannah.best.hits.seqs           <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Ophiophagus.hannah, input.seqs=Ophiophagus.hannah.genome_url, output.path="./Ophiophagus.hannah.UCEs.fas")
Pantherophis.guttatus.best.hits.seqs        <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Pantherophis.guttatus, input.seqs=Pantherophis.guttatus.genome_url, output.path="./Pantherophis.guttatus.UCEs.fas")
Protobothrops.mucrosquamatus.best.hits.seqs <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Protobothrops.mucrosquamatus, input.seqs=Protobothrops.mucrosquamatus.genome_url, output.path="./Protobothrops.mucrosquamatus.UCEs.fas")
Python.bivittatus.best.hits.seqs            <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Python.bivittatus, input.seqs=Python.bivittatus.genome_url, output.path="./Python.bivittatus.UCEs.fas")
Vipera.berus.best.hits.seqs                 <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Vipera.berus, input.seqs=Vipera.berus.genome_url, output.path="./Vipera.berus.UCEs.fas")
Thamnophis.sirtalis.best.hits.seqs          <- REEs::get.seqs.from.blastTable(input.blastTable=UCEs.best.hits.Thamnophis.sirtalis, input.seqs=Thamnophis.sirtalis.genome_url, output.path="./Thamnophis.sirtalis.UCEs.fas")
```

Output sequences in fasta format can be downloaded here: [Crotalus.horridus.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.horridus.UCEs.fas), [Crotalus.mitchellii.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.mitchellii.UCEs.fas), [Ophiophagus.hannah.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus.hannah.UCEs.fas), [Pantherophis.guttatus.UCEs.fasta](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis.guttatus.UCEs.fas), [Protobothrops.mucrosquamatus.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops.mucrosquamatus.UCEs.fas), [Python.bivittatus.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Python.bivittatus.UCEs.fas), [Thamnophis.sirtalis.UCEs.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis.sirtalis.UCEs.fas), [Vipera.berus.UCEs.fasta](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Vipera.berus.UCEs.fas). **Important note**: Sequence headers in the output fasta files (the text on lines starting with ">") have the format "querySequenceName_Subject=subjectContigName_ContigStart_ContigEnd". For example, the *Crotalus horridus* sequence that is the best match to UCE 1003 of *Micrurus fulvius* has the header "micrurus_fulvius_uce-1003_Subject=LVCR01005387.1_10818_10626". This is not intuitive because the only species in the header is *Micrurus fulvius*, even though the sequence is actually from *Crotalus horridus*.

To avoid confusion, I renamed sequences to the more intuitive header format "Genus_species_uce-XXXX", which matches the format used by Streicher and Weins (2017), where XXXX is a ID number for the UCE. For example: "Crotalus_horridus_uce-1003" is used as the name for the homolog of "micrurus_fulvius_uce-1003".

```
### Renaming each species' UCE sequences to match the format of sequence names used by Streicher and Wiens (2017): "Genus_species_uce-XXXX", where XXXX is a number.

names(Crotalus.horridus.best.hits.seqs)            <- mgsub(c("micrurus_fulvius_","_Subject.+"), c("Crotalus_horridus_",""), names(Crotalus.horridus.best.hits.seqs))
names(Crotalus.mitchellii.best.hits.seqs)          <- mgsub(c("micrurus_fulvius_","_Subject.+"), c("Crotalus_mitchellii_",""), names(Crotalus.mitchellii.best.hits.seqs))
names(Ophiophagus.hannah.best.hits.seqs)           <- mgsub(c("micrurus_fulvius_","_Subject.+"), c("Ophiophagus_hannah_",""), names(Ophiophagus.hannah.best.hits.seqs))
names(Pantherophis.guttatus.best.hits.seqs)        <- mgsub(c("micrurus_fulvius_","_Subject.+"), c("Pantherophis_guttatus_",""), names(Pantherophis.guttatus.best.hits.seqs))
names(Protobothrops.mucrosquamatus.best.hits.seqs) <- mgsub(c("micrurus_fulvius_","_Subject.+"), c("Protobothrops_mucrosquamatus_",""), names(Protobothrops.mucrosquamatus.best.hits.seqs))
names(Python.bivittatus.best.hits.seqs)            <- mgsub(c("micrurus_fulvius_","_Subject.+"),c("Python_bivittatus_",""), names(Python.bivittatus.best.hits.seqs))
names(Vipera.berus.best.hits.seqs)                 <- mgsub(c("micrurus_fulvius_","_Subject.+"), c("Vipera_berus_",""), names(Vipera.berus.best.hits.seqs))
names(Thamnophis.sirtalis.best.hits.seqs)          <- mgsub(c("micrurus_fulvius_","_Subject.+"), c("Thamnophis_sirtalis_",""), names(Thamnophis.sirtalis.best.hits.seqs))

### Save a copy of each of each species' renamed UCE sequences
writeXStringSet(Crotalus.horridus.best.hits.seqs,"Crotalus.horridus.UCEs_renamed.fas")
writeXStringSet(Crotalus.mitchellii.best.hits.seqs,"Crotalus.mitchellii.UCEs_renamed.fas")
writeXStringSet(Ophiophagus.hannah.best.hits.seqs,"Ophiophagus.hannah.UCEs_renamed.fas")
writeXStringSet(Pantherophis.guttatus.best.hits.seqs,"Pantherophis.guttatus.UCEs_renamed.fas")
writeXStringSet(Protobothrops.mucrosquamatus.best.hits.seqs,"Protobothrops.mucrosquamatus.UCEs_renamed.fas")
writeXStringSet(Python.bivittatus.best.hits.seqs,"Python.bivittatus.UCEs_renamed.fas")
writeXStringSet(Thamnophis.sirtalis.best.hits.seqs,"Thamnophis.sirtalis.UCEs_renamed.fas")
writeXStringSet(Vipera.berus.best.hits.seqs,"Vipera.berus.UCEs_renamed.fas")
```
The renamed sequences can be downloaded here: [Crotalus.horridus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.horridus.UCEs_renamed.fas), [Crotalus.mitchellii.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Crotalus.mitchellii.UCEs_renamed.fas), [Ophiophagus.hannah.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus.hannah.UCEs_renamed.fas), [Pantherophis.guttatus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis.guttatus.UCEs_renamed.fas), [Protobothrops.mucrosquamatus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops.mucrosquamatus.UCEs_renamed.fas), [Python.bivittatus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Python.bivittatus.UCEs_renamed.fas), [Thamnophis.sirtalis.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis.sirtalis.UCEs_renamed.fas), [Vipera.berus.UCEs_renamed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/UCEs.In.Snake.Genomes/Vipera.berus.UCEs_renamed.fas)


**5**. I used the function REEs::align.shared.loci to align the set of UCEs found in all snake genomes; this function invokes MAFFT for multisequence alignment.

```
### Define an character vector of paths to the input (unaligned) UCE sequence files. These are the sequences obtained from step 4. Note: Using the original or "renamed" (with headers renamed) sequences produces the same output.
UCEs.paths <- c("Crotalus.horridus.UCEs_renamed.fas","Crotalus.mitchellii.UCEs_renamed.fas","Ophiophagus.hannah.UCEs_renamed.fas","Pantherophis.guttatus.UCEs_renamed.fas","Protobothrops.mucrosquamatus.UCEs_renamed.fas","Python.bivittatus.UCEs_renamed.fas","Thamnophis.sirtalis.UCEs_renamed.fas","Vipera.berus.UCEs_renamed.fas")

### Align the UCEs that are found in each of the files specified by UCEs.paths
aligned.UCEs <- align.shared.loci(input.seqs=UCEs.paths, indv=c("Pantherophis.guttatus", "Thamnophis.sirtalis", "Vipera.berus", "Python.bivittatus", "Protobothrops.mucrosquamatus", "Ophiophagus.hannah", "Crotalus.mitchellii", "Crotalus.horridus"), reference.indv=7, seqname.str.delim="_", seqname.str.loc=3, output.dir="./MAFFT-aligned-UCEs")

```
Result: 2,968 of the UCEs from *Micrurus fulvius* were found in all of the snake genomes. The alignments for these shared UCEs can be download here: [MAFFT-aligned-UCEs.zip](https://osf.io/wpcr7/download).


**6**. Extract *Thamnophis sirtalis* sequences and filter by minimum size. Sort the remaining sequences by UCE name and then extract the middle 161nt of the first 1,000 loci.

```
### Character vector of paths to the UCE alignments
UCE.alignment.filenames    <- list.files(path="~/MAFFT-aligned-UCEs",full.names=T)

### Character vector of locus names (UCE ID number), with format "UCE.XXXX"
UCE.shortnames             <- gsub("uce-","UCE.",REEs::nameFromPath(UCE.alignment.filenames,include.extension=F),ignore.case=T)

### Read in UCE alignments (if they arent already loaded).
aligned.UCEs <- lapply(X=UCE.alignment.filenames,FUN=function(x){Biostrings::readDNAStringSet(x)})

### Set alignment names to values in UCE.shortnames
names(aligned.UCEs) <- UCE.shortnames

### Sort the list of alignments by alignment name
alignments.sorted <- aligned.UCEs[sort(UCE.shortnames)]

### Get Thamnophis sirtalis sequence from the filtered alignments
Thamnophis.sirtalis.seqs         <- REEs::collapse.DNAStringSet(lapply(X=alignments.sorted,FUN=function(X){X["Thamnophis.sirtalis"]}),use.seqnames=F)

### Remove gaps from Thamnophis sirtalis sequences
Thamnophis.sirtalis.noGaps.seqs <- DECIPHER::RemoveGaps(Thamnophis.sirtalis.seqs)

### Filter UCEs if T. sirtalis sequence length is less than 200nt
Thamnophis.sirtalis.UCEs.200 <- Thamnophis.sirtalis.noGaps.seqs[which(width(Thamnophis.sirtalis.noGaps.seqs)>=200)]

### Keep first 1000 UCEs of Thamnophis.sirtalis.UCEs.200
Thamnophis.sirtalis.UCEs.200.1000 <- Thamnophis.sirtalis.UCEs.200[1:1000]

### Get the middle 160nt (actually 161nt) from each sequence in Thamnophis.sirtalis.UCEs.200.1000. These are the UCE sequences submitted to Arbor for ultrastringent filtration and probe design. SnakeCap Probes are 120nt – choosing slightly larger targets allowed for some adjustments when designing probes.
# To get the position of the middle base of sequences in Thamnophis.sirtalis.UCEs.200.1000, divide sequence lengths by 2 (rounding down) 
middle.base   <- ceiling(width(Thamnophis.sirtalis.UCEs.200.1000)/2)
# Define the start and end positions of targets as 80nt upstream and downstream of the middle base.
start.target  <- middle.base-80
end.target    <- middle.base+80
# Define an IRanges and GRanges objects defining the UCE target ranges in Thamnophis.sirtalis.UCEs.200.1000
subranges     <- IRanges::IRanges(start=start.target ,end=end.target,names=names(Thamnophis.sirtalis.UCEs.200.1000))
gsubranges    <- GenomicRanges::GRanges(seqnames=names(Thamnophis.sirtalis.UCEs.200.1000),ranges=subranges)
# Extract UCE targets from Thamnophis.sirtalis.UCEs.200.1000
target.UCEs   <- Biostrings::getSeq(Thamnophis.sirtalis.UCEs.200.1000, gsubranges)

### Save UCE targets
Biostrings::writeXStringSet(x=target.UCEs,filepath="./targetUCEs.1000.fas",format="fasta")
```

Result: 1,000 *T. sirtalis* UCEs chosen as targets: [targetUCEs.1000.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/UCEs/targetUCEs.1000.fas). These correspond to WeinellEntry targets 2153–3152 in [Version2_additional-targets_Entry1899to3152_20Sep2018.fasta](https://github.com/JeffWeinell/SnakeCap/raw/main/ArborFiles/Version2_additional-targets_Entry1899to3152_20Sep2018.fasta). After ultrastringent filtering, probes were designed for 907 of the 1,000 UCEs submitted to Arbor Biosciences. See [ultra-stringent filtering](#ultrastringentFiltering) section.


<a name="Methods.SelectingddRAD"></a>
## Selecting the set of target ddRAD-like loci

#### Overview:

Identify and obtain all regions of the *Thermophis baileyi* genome that begin with the recognition sequence of the restriction enzyme (RE) *Sbfi* and end with the recognition sequence of the RE *EcoR1*, and which have a length of 900-1000nt. Then I used BLASTn to search for homologs of the *Thermophis baileyi* ddRAD-like loci in the *T. sirtalis* genome. I filtered the set of ddRAD-like loci to include only those with a single, high quality match found in *T. sirtalis* genome. I trimmed the ends of sequences if there was a string of ambiguous bases nearby. A subset of the retained loci were submitted to Arbor for ultrastringent filtering and then probe design.

<!---
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
--->

#### Details:

Step-by-step description of pipeline used to select ddRAD-like loci to target:

**1**. I used the function REEs::proposeLoci.ddRADlike to propose a set of target loci ("ddRAD-like" loci) that are expected to be distributed relatively randomly throughout the genome, and which have lengths between 900-1000nt. These loci ddRAD-like because the method involves searching for the recognition sequences of two restriction enzymes. The function returns a table containing the genomic coordinates of all regions that begin with the recognition sequence of a restriction enzyme A (= SbfI: "CCTGCAGG") and end with the recognition sequence of a restriction enzyme B (= EcoR1: "GAATTC"), and which have a sequence length that is within a specified range (900-1000nt). Important Note: the function proposeLoci.ddRADlike should not be thought of as a way to conduct in silico ddRAD, and it is possible that the proposed loci contains internal RE recognition sites, which wouldn't happen in real ddRAD.

```
### Define the recognition sequence for SbfI restriction enzyme
SbfI.Seq                              <- "CCTGCAGG"
### Define the recognition sequence for EcoR1 restriction enzyme
EcoR1.Seq                             <- "GAATTC"
### Define URL path to Thermophis baileyi genome
Thermophis.baileyi_genome.url         <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/457/575/GCA_003457575.1_DSBC_Tbai_1.0/GCA_003457575.1_DSBC_Tbai_1.0_genomic.fna.gz"

### Use the function proposeLoci.ddRADlike (REEs package) to find 900-1000nt regions of contigs that begin with one of the two recognition sequences (SbfI.Seq or EcoR1.Seq) and end with the other recognition sequence.
proposed.ddRAD.loci.coordinates       <- REEs::proposeLoci.ddRADlike(input.seqs=Thermophis.baileyi_genome.url,output.dir="/ddRAD-like/",recognitionSeqs=c(SbfI.Seq,EcoR1.Seq),lim.lengths=c(900,1000),save.tables=T)
```
Result: Table containing coordinates for 2,337 proposed ddRAD-like loci: [ProposedLoci_CCTGCAGG-GAATTC_output_900to1000.txt](https://github.com/JeffWeinell/SnakeCap/raw/main/ddRAD/ProposedLoci_CCTGCAGG-GAATTC_output_900to1000.txt).

**2**. I used the function get_ncbi_sequences (REEs package) to obtain the sequences of the loci proposed in step 1. For loci on the antisense strands of contigs I first downloaded the corresponding reverse complement sequence (from the sense strand), and then I used the function reverseComplement (Biostrings package) to obtain the correct antisense strand target.
```
### Defining vectors containing the Genbank accession and range of positions to download. This information is from the second, third, and fourth columns of the table that was generated by the function proposeLoci.ddRADlike. Note that the start and end positions are relative to the sense strand of the contig. The coordinates of proposed loci on the antisense strand have the format "ContigAccession:cEndPosition-StartPosition", which is equal to the reverse complement of a sequence with coordinates "ContigAccession:StartPosition-EndPosition".
ContigAccession <- proposed.ddRAD.loci.coordinates$ContigAccession
LocusStart      <- proposed.ddRAD.loci.coordinates$StartPosition
LocusEnd        <- proposed.ddRAD.loci.coordinates$EndPosition

### Download genome to temporary file, and then extract and import from positions LocusStart to LocusEnd of the sense strand of the respective contig in ContigAccession. The "outfile" parameter is set to NULL because we need to take the reverse complement of some of the sequences before saving. 
Thermophis.ddradlike.proposed.seqs <- REEs::get_ncbi_sequences(outfile=NULL,input.seqs=Thermophis.baileyi_genome.url,accessionList=ContigAccession,startList=LocusStart,endList=LocusEnd,strandList="1",db="nuccore",rettype="fasta",retmode="text",if.outside.range="partial",trim.ambiguous = FALSE)

### Numerical vector indicating which of the proposed loci are on the antisense strand. These are the loci that need to be reverse complemented.
which.reverse.complement <- grep(":c",proposed.ddRAD.loci.coordinates$coordinates)

### Reverse complement the "ith" sequence in Thermophis.ddradlike.seqs for each i in the vector which.reverse.complement
Thermophis.ddradlike.seqs[which.reverse.complement] <- Biostrings::reverseComplement(Thermophis.ddradlike.seqs[which.reverse.complement])

### Update sequence names to account for reverse complementing. The new sequence names are the character strings in the coordinates column (first column) of the table proposed.ddRAD.loci.coordinates.
names(Thermophis.ddradlike.seqs) <- proposed.ddRAD.loci.coordinates$coordinates

### Save sequences.
writeXStringSet(Thermophis.ddradlike.seqs,"Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000.fas")
```
*Thermophis baileyi* sequences for proposed ddRAD-like loci can be downloaded here: [Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/ddRAD/Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000.fas). **Note: QLTV01020412.1:1-995, which was targetted, is not in the table of proposed loci (nor is any locus on that contig); need to figure out why.**

**3**. Trim loci that have strings of ambiguous bases to keep the longest non-ambiguous region.
```
### If not already loaded, read in the set of proposed loci sequences from step two.
proposed.loci <- Biostrings::readDNAStringSet("Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000.fas")

### Vector indicating which loci contain a string of ambiguous bases (Ns).
which.Ns <- grep("N+N",proposed.loci)

### DNAStringSet object containing the loci to trim
loci.to.trim <- proposed.loci[which.Ns]

### Start and end positions of the strings of Ns
N.locs <- stringr::str_locate(string=loci.to.trim,pattern="N+N")

### Start and end positions of the 5' and 3' non-ambiguous strings.
left.pos     <- cbind(1,(N.locs[,1]-1))
right.pos    <- cbind((N.locs[,2]+1),width(loci.to.trim))

### Length of 3' and 5' non-ambiguous regions
length.left  <- (N.locs[,1]-1)
length.right <- width(loci.to.trim)-N.locs[,2]
which.longer <- apply(X=cbind(length.left,length.right),MARGIN=1,FUN=function(x){which(x==max(x))})
which.left   <- which(which.longer==1)
which.right  <- which(which.longer==2)

### DNAStringSets containing the trimmed loci.
loci.trimmed.left  <- subseq(loci.to.trim[which.left],start=rep(1,length(which.left)),end=left.pos[,2][which.left])
loci.trimmed.right <- subseq(loci.to.trim[which.right],start=right.pos[,1][which.right],end=right.pos[,2][which.right])

### Old contig coordinates of loci to be trimmed
old.coordinates.left        <- mat.strsplit(mgsub(c(".+:c",".+:"),c("",""),names(loci.trimmed.left)),split="-")
old.coordinates.right       <- mat.strsplit(mgsub(c(".+:c",".+:"),c("",""),names(loci.trimmed.right)),split="-")
mode(old.coordinates.left)  <- "numeric"
mode(old.coordinates.right) <- "numeric"

### New contig coordinates of loci to be trimmed
new.coordinates.left        <- cbind(old.coordinates.left[,1],old.coordinates.left[,1]+(left.pos[which.left,2]-1))
new.coordinates.left[grep(":c",names(loci.trimmed.left)),2] <- old.coordinates.left[grep(":c",(names(loci.trimmed.left))),1]-(left.pos[which.left[grep(":c",(names(loci.trimmed.left)))],2]-1)

new.coordinates.right                                         <- cbind((old.coordinates.right[,1]+right.pos[which.right,1]-1),old.coordinates.right[,2])
new.coordinates.right[grep(":c",names(loci.trimmed.right)),1] <- old.coordinates.right[grep(":c",(names(loci.trimmed.right))),1]-(right.pos[which.right[grep(":c",(names(loci.trimmed.right)))],1]-1)

new.names.left        <- paste0(mgsub(c(":c[1-9].+",":[1-9].+"),c(":c",":"),names(loci.trimmed.left)),new.coordinates.left[,1],"-",new.coordinates.left[,2])
new.names.right       <- paste0(mgsub(c(":c[1-9].+",":[1-9].+"),c(":c",":"),names(loci.trimmed.right)),new.coordinates.right[,1],"-",new.coordinates.right[,2])

### Rename the trimmed loci to the new genomic coordinates
names(loci.trimmed.left)  <- new.names.left
names(loci.trimmed.right) <- new.names.right
loci.trimmed <- c(loci.trimmed.left,loci.trimmed.right)
loci.trimmed <- loci.trimmed[match(unique(names(loci.trimmed)),names(loci.trimmed))]

### DNAStringSet object of the proposed loci, including the trimmed loci
proposed.loci.trimmed                        <- c(proposed.loci[-which.Ns],loci.trimmed)

### Save proposed.loci.trimmed
writeXStringSet(proposed.loci.trimmed,"Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000_trimmed.fas")

```
Proposed ddRAD-like loci after trimming can be downloaded here: [Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000_trimmed.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/ddRAD/Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000_trimmed.fas).

**4**. I randomly selected loci from the set of trimmed, proposed ddRAD-like loci until all remaining baits (n = ) in the 20K bait kit were assigned. **Finish adding the code used to calculate baits per locus**
```
### Read in the trimmed loci if not already loaded
proposed.loci.trimmed             <- Biostrings::readDNAStringSet("Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000_trimmed.fas")

### Randomly shuffle the order of the loci
proposed.loci.trimmed.randomOrder <- proposed.loci.trimmed[sample(1:length(proposed.loci.trimmed),n=lenth(proposed.loci.trimmed),replace=F)]

### Calculate the number of baits needed to target each locus in proposed.loci.trimmed.randomOrder
baits.per.locus <- 

### Create a vector of the rolling sum of baits per locus
baits.per.locus.rollingSum <- REEs::rollSum(baits.per.locus)

### Keep the first N loci in proposed.loci.trimmed.randomOrder in which rolling sum of baits per locus is closest to but not exceeding the number of unassigned baits remaining in the bait kit.

ddrad.targets.seqs <- proposed.loci.trimmed.randomOrder[which(baits.per.locus.rollingSum < 3,695)]

### Save the set of proposed ddRAD-like loci to target
writeXStringSet(ddrad.targets.seqs,"/Users/alyssaleinweber/Documents/SequenceCapture-GitHub/ddRAD-like/Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000_trimmed_targetted.fas")
```
Set of ddRAD-like loci to target can be downloaded here: [Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000_trimmed_targetted.fas](https://github.com/JeffWeinell/SnakeCap/raw/main/ddRAD/Thermophis_ProposedLoci_CCTGCAGG-GAATTC_900to1000_trimmed_targetted.fas).
**Note: QLTV01020412.1:1-995, which was also targetted, is not in this file**


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

**Table 6**. Summary of the MHC filtered hit table. The table below shows the number MHC loci (columns 2-5) with at least n matches (column 1) in the genome of *T. sirtalis* with at least *x* percent coverage and percent identical sites across the covered region.
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
 
The remaining 27 MHC loci (non-REEs) and six others that were targetted as REEs (entries 248, 559, 728, 787, 891, and 1944) were included in the final set of target loci for which baits were synthesized. WeinellEntry1944 was later recognized as a duplicate of the immune target WeinellEntry1857.

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
### Selecting vision loci:

I used blastn to search for the vision loci probes from Schott et al. (2017) (which were from *Anolis*, *Columba*, *Gallus*, and *Pelodiscus*, *Sceloporus*, or *Python*) within the snake genomes. Most of the SnakeCap probes for these loci are designed from *Ophiophagus* (n = 88), but some probes were designed from *Thamnophis* (n = 21), *Protobothrops* (n = 5), *Pantherophis* (n = 3), or *Python* (n = 2), when blastn of Schott et al 2017 probes did not yield a strong match in *Ophiophagus*.

<a name="ultrastringentFiltering"></a>
<a name="ProbeSynthesis"></a>
## Ultra-stringent filtration and probe synthesis

After choosing the target loci, probes were designed by Arbor Biosciences with the following specifications: 50% tiling, 120nt/probe; 20,020 probes in total. See **Target-loci_Coverage_graph_22October2020.pdf** for a visual summary of target loci, probes, probe coverage, and features of loci including genes, mRNA/transcribed regions, and protein-coding (CDS) regions. This graph was generated with **graph_target_and_features.R** and then filesize reduction in Adobe Acrobat.

**Table 7.** Genomes from which synthesized baits were designed from.
species | genome assembly accession | contig prefixes | n target loci with baits designed from species | n bait-containing contigs
----|----|----|----|----
*Thamnophis sirtalis* | GCF_001077635.1 | NW_0136; LFLD01 | 2,619  | 1,412 
*Ophiophagus hannah* | GCA_000516915.1  | AZIM01 | 151 | 86 
*Crotalus horridus* | GCA_001625485.1 | LVCR01 | 2  | 2 
*Python bivittatus* | GCF_000186305.1 | NW_0065; AEQU02 | 18  | 7 
*Probothrops mucrosquamatus* | GCF_001527695.1 | NW_0153; BCNE02 | 9  | 9
*Pantherophis guttatus* | GCA_001185365.1 | JTLQ01 | 2  | 2
*Thermophis baileyi* | GCA_003457575.1 | QLTV01 | 328  | 229

<a name="References"></a>
## References:

Bankevich A., et al. 2012. Spades: a new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology* 19, 455–477. https://doi.org/10.1089/cmb.2012.0021

Armstrong J., Hickey G., Diekhans M. et al. 2020. Progressive Cactus is a multiple-genome aligner for the thousand-genome era. *Nature* 587, 246–251. https://doi.org/10.1038/s41586-020-2871-y

Bushnell B. 2014 Bbmap: a fast, accurate, splice-aware aligner. Lawrence Berkeley National Lab (LBNL). https://www.osti.gov/biblio/1241166

Bushnell B., Rood J., and Singer E. 2017. Bbmerge: accurate paired shotgun read merging via overlap. PLoS One 12, e0185056. https://doi.org/10.1371/journal.pone.0185056

Hickey G., Paten B., Earl D., Zerbino D., and D. Haussler. 2013. HAL: a hierarchical format for storing and analyzing multiple genome alignments, Bioinformatics, 29(10) 1341–1342. https://doi.org/10.1093/bioinformatics/btt128

Holthaus K.B., Mlitz V., Strasser B., Tschachler E., Alibardi L., and L. Eckhart. 2017. Identification and comparative analysis of the epidermal differentiation complex in snakes. *Scientific Reports* 7, 45338. https://doi.org/10.1038/srep45338

Hutter C.R., Cobb K.A., Portik D., Travers S., Wood Jr. P.L., and R.M. Brown. 2019. FrogCap: A modular sequence capture probe set for phylogenomics and population genetics for Anurans, assessed across multiple phylogenetic scales. *bioRxiv* 825307. https://doi.org/10.1101/825307

Schott R.K., Panesar B., Card D.C., Preston M., Castoe T.A., and B.S.W. Chang. 2017. Targeted Capture of Complete Coding Regions across Divergent Species. *Genome Biology and Evolution* 9(2), 398–414. http://doi.org/10.1093/gbe/evx005

Streicher J.W., and J.J Wiens J.J. 2017. Phylogenomic analyses of more than 4000 nuclear loci resolve the origin of snakes among lizard families. *Biology Letters* 13, 20170393. http://doi.org/10.1098/rsbl.2017.0393

