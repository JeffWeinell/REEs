# SnakeCap Sequence Capture Probe Set:

## Contents

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
## Description of SnakeCap Probe Set

The probe set includes 20,020 probes for 3,129 single-copy loci (1,517,011 nt) shared across snakes. The target loci are categorized into four types: (1) rapidly evolving exons (REEs; n = 1,653), (2) ultra-conserved elements (UCEs; n = 907), (3) ddRAD-like loci (n = 328), (4) and functionally interesting genes, which includes 27 major histocompatibility complex (MHC) genes, 119 vision genes, and 95 scalation genes.

REEs include one or more entire exons and one or both exon-flanking regions, and range in length from 121 to 7,501 nt. I used a modified version of the FrogCap pipeline (Hutter et al., 2019) to select the optimal set of REEs from an alignment of snake exomes.

SnakeCap UCEs are a subset of the *Micrurus fulvius* UCEs from Streicher and Wiens (2017).

ddRAD-like loci are shared, single-copy loci identified from in-silico ddRAD using recognition sites for SbfI and EcoRI restriction enzymes.

Functional loci included entire or partial gene regions that have previously been predicted or known to function in either (1) vertebrate immune systems, (2) vision, (3) or scalation.

#### Table 2.1 For each type of locus: genomic region targeted, number of loci (nloci), total number of nucleotides targeted (nnt), and nucleotide lengths (nt/locus) of the shortest and longest loci.

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
## Methods

<a name="Methods.SelectingREEs"></a>
### Selecting the set of target REEs

<!-- <a name="Methods.SelectingREEs.overview"></a> -->
#### Overview: 

<!-- <a name="Methods.SelectingREEs.detailed"></a> -->
#### Detailed, step-by-step methods for how I chose the set of target REEs

1. I downloaded the [*Thamnophis sirtalis* genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz) and its associated annotation table: [GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz) (n = 559,130 features annotated). Then, I renamed the contigs in the genome file to have the following format: **Thamnophis_sirtalis_GCF_001077635.1_read1**, **Thamnophis_sirtalis_GCF_001077635.1_read2**, etc., and saved this renamed genome in sequential fasta format: [ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3.zip](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/exomes/ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3.zip). The two-column, tab-delimited table [Scaffold-Name-Key.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/main/exomes/Scaffold-Name-Key.txt?token=AJJOG2UQ6MDA7UY2U4R6BFS7ZDZYS) includes the new contig name in the first column and the original contig name in the second column:

```
### URL to the Thamnophis sirtalis genome feature table.
Thamnophis.sirtalis_GFF.url       <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"

### URL to the Thamnophis sirtalis genome (fasta formatted sequences).
Thamnophis.sirtalis_genome.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz"

### Loading the feature table into R
Thamnophis.sirtalis_GFF           <- data.table::fread(Thamnophis.sirtalis_GFF.url,col.names=c("seqname","source","feature","start","end","score","strand","frame","attribute"),skip=8,sep="\t",fill=TRUE,blank.lines.skip=TRUE)

### Removing incomplete/empty feature lines:
Thamnophis.sirtalis_GFF <- Thamnophis.sirtalis_GFF[-which(Thamnophis.sirtalis_GFF$feature=="")]

### Saving the feature table as .tsv format (i.e., without the GFF headers or incomplete features):
write.table("~/ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3",sep="\t")
```

<!---
At some point this file becomes involved: ref_Thamnophis_sirtalis-6.0_top_level.gff3
--->
<!---
ScaffoldName | RefSeq.ScaffoldAccession
------------ | -------------
Thamnophis_sirtalis_GCF_001077635.1_read1 | RefSeq accession number (the ID from column 1 of the annotation table)
Thamnophis_sirtalis_GCF_001077635.1_read2 | ...
...|...
Thamnophis_sirtalis_GCF_001077635.1_readN|...
--->

2. To extract the *Thamnophis sirtalis* exome from the genome, I first used the function **filter.annotationTable** to: (1) filter the original annotation table to only include annotations for regions that are both CDS regions and ≥ 120bp in length (result: n = 115,907 annotations for 77,329 unique features), and (2) write the filtered annotation table to a file: **CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_withGenBankAcc_longer120bp.gff3**. Then, I used the function **get.exome.from.annotationTable** which extracts the DNA sequences included in the filtered annotation table from the genome. The extracted DNA sequences (the exome) were written in fasta format to: **Thamnophis_sirtalis_exome_longer120bp.fas**.

```
### Filtering the feature table to only include CDS features at least 120bp.
Thamnophis.sirtalis_GFF_CDS_longer120bp <- filter.annotationTable(input.gff=Thamnophis.sirtalis_GFF,output.gff="CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_longer120bp.gff3",feature.type="CDS",min.length=120)

### Extracting sequences in the feature table Thamnophis.sirtalis_GFF_CDS_longer120bp from the genome.
get.exome.from.annotationTable(species.name="Thamnophis_sirtalis",genome.filepath="~/Thamnophis_sirtalis_GCF_001077635.1_genome_renamed_sequential.fas",input.gff="CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_longer120bp.gff3",output.dir.exome = "~/exomes/",additional.ID="Scaffold-Name-Key.txt")
```

3. I downloaded all squamate genomes available from NCBI and then searched within each for exons homologous to those in the *T. sirtalis* exome (**Table 4**).

**Table 4**. Genomes used to select REEs included all squamate genomes available from NCBI.
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

I queried each *T. sirtalis* exon against each squamate genome using **tblastx** (allowing up to 50 matches to be saved per query). This step produces a Hit Table for each species with stats on each query/target match, including the bitscore, which is a measurement how good the match is (how likely the match corresponds to homology). **Note to self**: I used the cluster submission file **TBLASTX.sh** to perform this step.

Then, I filtered each of the full (i.e., 50 matches/query) hit tables to include only the best match/query (max bitscore) using the R function **reportBestMatches**.

```
reportBestMatches(input.dir="TBlastXResults/",species.names=c("Anolis_carolinensis","Crotalus_mitchellii","Pogona_vitticeps","Gekko_japonicus","Ophiophagus_hannah","Vipera_berus","Crotalus_horridus","Thamnophis_sirtalis","Python_bivittatus","Protobothrops_mucrosquamatus","Pantherophis_guttatus"),output.dir=NA,blastMethod="tblastx",locusType="exons")
```

4. Extract exomes (minimum exon length 120nt) for each squamate species. For each species and exon, I extracted the DNA sequence of the best match in the filtered hit table from step 2, and I saved these sequences to a fasta file; this step was performed using the function **get.exome.from.blastTable**.

```
### Example for Vipera berus
get.exome.from.blastTable(species="Vipera_berus",genome.filepath="~/Vipera_berus_GCA_000800605.1_Vber.be_1.0_genomic.fna",input.blastTable="Vipera_berus.tblastx.exons.best.txt",output.dir.exome="~/exomes/")
```

5. Align shared exons and calculate stats. I used the R function **makeExomeStatsTable** to do all of the following:
  - obtain the set of *T. sirtalis* exons present in all exomes (from step 3), i.e., the shared exons
  - perform multiple sequence alignment (MAFFT algorithm) for each of the shared exons
  - calculate a set of stats for each shared exon alignment, and save results to the file **stats_exome_data_TBLASTX.txt** (results: stats for 66,489 alignments). 
  
```
### "species" parameter includes the set of the squamates
### "subgroup" parameter includes the set of the snakes
### "is.primary.exome" was set to 1 because Thamnophis sirtalis is the first species in the "species" parameter vector, and this is the species that baits were designed from

stats_exome <- makeExomeStatsTable(exomes.filepaths=list.files(path="~/exomes/",full.names=T), annotationTable.path="CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_withGenBankAcc_longer120bp.gff3", species=c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus","Protobothrops_mucrosquamatus","Pantherophis_guttatus","Anolis_carolinensis","Pogona_vitticeps","Gekko_japonicus"), subgroup=c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus","Protobothrops_mucrosquamatus","Pantherophis_guttatus"), output.dir="~/exomes/", is.primary.exome=1)
```

The file **stats_exome_data_TBLASTX.txt** contains a table with NCBI annotation information for each exon. Rows correspond to exons, and columns include the following: (1) *Thamnophis sirtalis* NCBI Reference Sequence ID and location of exon (stop_start) **(start_stop?)** (2) number of species in alignment (always 11, because 11 species included, and only shared loci were aligned), (3) number of sites with at least four species represented, (4) number of parsimony informative sites, (5) percent of sites parsimony informative, (6) mean pairwise percent genetic similarity to *T. sirtalis*, (7–17) percent genetic similarity of each species to *T. sirtalis*, (18) exon alignment width, (19) width of *T. sirtalis* exon, (20) gene name associated with the exon (extracted from the last column of the filtered annotation table of *T. sirtalis*), (21) mean number of variable sites compared to *T. sirtalis*, (22) minimum percent genetic similarity to *T. sirtalis* (among the 11 species included), (23) minimum percent genetic similarity to *T. sirtalis* among the snakes included.

6. Filter loci to include the most rapidly evolving loci while considering constraints imposed by the bait kit. To do this, I used the function **pick.loci**, which: (1) filtered out loci if minimum percent genetic similarity (among snakes) to *T. sirtalis* was < 65% or = 100% (result: 64,546 exons pass step 1). (2) For genes with multiple exons, kept only the exon with the lowest mean pairwise genetic distance to *T. sirtalis* (16,650 exons pass step 2 and are saved to [stats_data_FastestExonPerGene_best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/REEs/stats_data_FastestExonPerGene_best.txt)). (3) Identified the set of exons maximizing the total number of variable sites (all exons in set) while meeting the following constraints and conditions: total number of nucleotides targetted ≤ 1.2Mb and number of baits ≤ 20,000 (these are constraints of the 20K my-baits kit), bait size = 120nt, and bait tiling = 50% overlap.

```
result.pident <- pick.loci(statsTable.path = "~/AlignedExonStats/stats_exome_data_TBLASTX.txt", output.dir = "~/AlignedExonStats/", primary.species = "Thamnophis_sirtalis", use.min.pident.subgroup = T, species.subgroup = c(7:14), min.pident.keep = c(65,100), max.capture.coverage = 1200000, write.stats.tables = T, plot.results = T, fast.stat = "pident")
```

Result = 2,068 REEs retained; the stats table for these loci was written to the file [stats_data_FastestExonPerGene_best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/REEs/stats_data_FastestExonPerGene_best.txt); an updated version of this stats table that includes the WeinellEntry locus names is [stats_data_FastestExonPerGene_best_20Nov2020.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/REEs/stats_data_FastestExonPerGene_best_20Nov2020.tsv).

These 2,068 REEs were submitted to Arbor Biosciences for probe design. Arbor performed ultrastringent filtration these loci which resulted in the removal of 70 REEs [Version1-loci-removed_ZeroBaitCoverageLoci.tsv](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/ArborFiles/Version1-loci-removed_ZeroBaitCoverageLoci.tsv), whereas 1,998 passed this step. An additional 123 REEs were removed because the baits designed to target these loci were all non-specific within the genomes of *T. sirtalis* and/or *Thermophis baileyi*; 212 other REEs were removed to allow some baits to be used to target other types of loci (UCEs, immune, scalation, vision, and ddRAD-like loci). The remaining 1,653 REEs were synthesized and included in the mybaits 20K bait kit (product no. 3001160). 

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
### Selecting the set of target UCEs

<!-- <a name="Methods.SelectingUCEs.overview"></a> -->
#### Overview:

Target UCEs include 907 of the 3,260 UCEs previously identified in *Micrurus fulvius* (Streicher and Wiens, 2017; **Table X**). First, I filtered the full set of *Micrurus* UCEs to only include those present in all NCBI snake genomes (n = 2,968 UCEs). Then, I filtered the shared set of UCEs to only include those with an alignment width > 200nt (2,551 UCEs retained; each UCE alignment included the eight snakes with published genomes). Next, I removed the following UCEs: uce-1843, uce-2179, uce-2433, uce-2465, uce-2498, uce-2890, uce-2960, and uce-3354 (not sure why I did this yet). I sorted the remaining 2,543 UCEs by UCE name and retained the first 1,000 loci in this set. Arbor Biosciences was able to synthesize probes for 907 of the 1,000 proposed target UCEs.

<!-- <a name="Methods.SelectingUCEs.detailed"></a> -->
#### Detailed, step-by-step methods for how I chose the set of target UCEs:

1. I downloaded the set of *Micurus fulvius* UCEs (n = 3,260) identified by Streicher and Wiens (2017). These were available as a fasta file called [micrurus_UCEs.fa](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/micrurus_UCEs.fa).

2. To find *Micrurus* UCEs in each of the other snake genomes, I queried each *Micrurus* UCE against each snake genome using blastn algorithm (saving ≤ 50 matches per query), which was implemented using the R wrapper function **blastnR**. Results were saved to the files: **Crotalus_horridus.blastn.UCEs.50hits.txt**, **Crotalus_mitchellii.blastn.UCEs.50hits.txt**, **Ophiophagus_hannah.blastn.UCEs.50hits.txt**, **Pantherophis_guttatus.blastn.UCEs.50hits.txt*, **Protobothrops_mucrosquamatus.blastn.UCEs.50hits.txt**, **Python_bivittatus.blastn.UCEs.50hits.txt**, **Thamnophis_sirtalis.blastn.UCEs.50hits.txt**, and **Vipera_berus.blastn.UCEs.50hits.txt**.

```
blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/GCA_000737285.1_CrotMitch1.0_genomic.fna",query.path="~/micrurus_UCEs.fa",output.path="~/Crotalus_mitchellii.blastn.UCEs.50hits.txt")

blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/Protobothrops_mucrosquamatus_GCF_001527695.2_P.Mucros_1.0_genomic.fna", query.path="~/micrurus_UCEs.fa",output.path="~/Protobothrops_mucrosquamatus.blastn.UCEs.50hits.txt")

blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/Pantherophus_guttatus_GCA_001185365.1_PanGut1.0_genomic.fna", query.path="~/micrurus_UCEs.fa",output.path="~/Pantherophis_guttatus.blastn.UCEs.50hits.txt")

blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/Python_bivittatus_GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna", query.path="~/micrurus_UCEs.fa",output.path="~/Python_bivittatus.blastn.UCEs.50hits.txt")

blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/Thamnophis_sirtalis_GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna" ,query.path="~/micrurus_UCEs.fa",output.path="~/Thamnophis_sirtalis.blastn.UCEs.50hits.txt")

blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/Ophiophagus_hannah_GCA_000516915.1_OphHan1.0_genomic.fna", query.path="~/micrurus_UCEs.fa",output.path="~/Ophiophagus_hannah.blastn.UCEs.50hits.txt")

blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/Vipera_berus_GCA_000800605.1_Vber.be_1.0_genomic.fna", query.path="~/micrurus_UCEs.fa",output.path="~/Vipera_berus.blastn.UCEs.50hits.txt")

blastnR(blastn.path="~/ncbi-blast-2.5.0+/bin/blastn",subject.path="~/Crotalus_horridus_GCA_001625485.1_ASM162548v1_genomic.fna", query.path="~/micrurus_UCEs.fa",output.path="~/Crotalus_horridus.blastn.UCEs.50hits.txt")
```

3. Then, I filtered the hit tables to include only the best match/query (max bitscore) using the R function **reportBestMatches**. Results were saved to the files: [Crotalus_horridus.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Crotalus_horridus.blastn.UCEs.best.txt), [Crotalus_mitchellii.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Crotalus_mitchellii.blastn.UCEs.best.txt), [Ophiophagus_hannah.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Ophiophagus_hannah.blastn.UCEs.best.txt) , [Pantherophis_guttatus.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Pantherophis_guttatus.blastn.UCEs.best.txt), [Protobothrops_mucrosquamatus.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Protobothrops_mucrosquamatus.blastn.UCEs.best.txt), [Python_bivittatus.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Python_bivittatus.blastn.UCEs.best.txt), [Thamnophis_sirtalis.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Thamnophis_sirtalis.blastn.UCEs.best.txt), [Vipera_berus.blastn.UCEs.best.txt](https://raw.githubusercontent.com/JeffWeinell/SnakeCap/blob/main/UCEs/UCEs.In.Snake.Genomes/Vipera_berus.blastn.UCEs.best.txt).

```
reportBestMatches(input.dir="~/BLAST_Micrurus-UCEs_vs_SnakeGenomes/",species.names=c("Crotalus_mitchellii","Crotalus_horridus","Ophiophagus_hannah","Pantherophis_guttatus","Protobothrops_mucrosquamatus","Python_bivittatus","Thamnophis_sirtalis","Vipera_berus"),output.dir=input.dir,blastMethod="blastn",locusType="UCEs")
```

4. To extract and save the set of best-match UCEs from each genome I used the function **get.UCEs.from.blastTable**.

```
get.UCEs.from.blastTable(species="Crotalus_horridus",genome.filepath=,input.blastTable="Crotalus_horridus.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")

get.UCEs.from.blastTable(species="Crotalus_mitchellii",genome.filepath=,input.blastTable="Crotalus_mitchellii.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")

get.UCEs.from.blastTable(species="Ophiophagus_hannah",genome.filepath=,input.blastTable="Ophiophagus_hannah.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")

get.UCEs.from.blastTable(species="Pantherophis_guttatus",genome.filepath=,input.blastTable="Pantherophis_guttatus.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")

get.UCEs.from.blastTable(species="Protobothrops_mucrosquamatus",genome.filepath=,input.blastTable="Protobothrops_mucrosquamatus.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")

get.UCEs.from.blastTable(species="Python_bivittatus",genome.filepath=,input.blastTable="Python_bivittatus.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")

get.UCEs.from.blastTable(species="Thamnophis_sirtalis",genome.filepath=,input.blastTable="Thamnophis_sirtalis.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")

get.UCEs.from.blastTable(species="Vipera_berus",genome.filepath=,input.blastTable="Vipera_berus.blastn.UCEs.best.txt",output.dir="~/UCEs.In.Snake.Genomes/")
```

5. To identify and align the set of UCEs found in all snake genomes, I used the function **align.bestHit.UCEs**. This function invokes MAFFT to perform multisequence alignment.

```
align.bestHit.UCEs(species.UCEs.filepaths=list.files(path="~/UCEs.In.Snake.Genomes/",full.names=T), output.dir="~/MAFFT-aligned-UCEs", species=c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus", "Protobothrops_mucrosquamatus","Pantherophis_guttatus"))
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
### Selecting the set of target ddRAD-like loci

#### Overview:

#### Detailed, step-by-step methods for how I chose the set of target ddRAD-like loci:

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
### Selecting MHC loci:

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
### Selecting scalation loci:

I targeted a subset of the genes included in the study by Holthaus et al. (2017). In that study, the authors identified homologous genes of the Epidermal Differentiation Complex (which are putatively involved in scalation) of *Python bivittatus* and *Ophiophagus hannah*. I downloaded the *Ophiophagus* scalation gene sequences using the table of genomic coordinates provided by Holthaus et al. (2017), and then used tblastn to search for and obtain homologous loci in *T. sirtalis*, *Protobothrops mucrosquamatus*, and *Crotalus horridus*.


<a name="Methods.SelectingVision"></a>
### Selecting vision loci:

I used blastn to search for the vision loci probes from Schott et al. (2017) (which were from *Anolis*, *Columba*, *Gallus*, and *Pelodiscus*, *Sceloporus*, or *Python*) within the snake genomes. Most of the SnakeCap probes for these loci are designed from *Ophiophagus* (n = 88), but some probes were designed from *Thamnophis* (n = 21), *Protobothrops* (n = 5), *Pantherophis* (n = 3), or *Python* (n = 2), when blastn of Schott et al 2017 probes did not yield a strong match in *Ophiophagus*.


<a name="ProbeSynthesis"></a>
### Probe Synthesis

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
### Taxa sampled

The following species were sequenced: *Achalinus spinalis*, *Aparallactus capensis*, *Aplopeltura boa*, *Elapsoidea sundevalli*, *Boiga irregularis*, *Buhoma depressiceps*, *Cerberus schneideri*, *Chamaelycus fasciatus*, *Cyclocorus lineatus*, *Oxyrhabdium* cf. *modestum*, *Oxyrhabdium leporinum*, *Oxyrhabdium modestum*, *Prosymna visseri*, *Psammodynastes pulverulentus*, *Pseudaspis cana*, *Pseudoxyrhopus tritaeniatus*, *Rhamphiophis oxyrhynchus*, *Scolecophis atrocinctus*, *Tantilla taeniata*.

<a name="LibraryPrep"></a>
### Sequence capture library prep

Conducted by Arbor Biosciences; eight samples/pool; ...

<a name="DNASequencing"></a>
### DNA sequencing

Novogene Illumina HiSeqX; paired-end sequencing, read length 150nt, insert size 400nt?.

<a name="PostSequencing"></a>
### Post-sequencing

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
- **01_Pre_Process_Reads_Apr10.R**
- **02_Assemble_Spades_Apr18.R**
- **03_Target-loci_matching_20Feb2020.R** (= **03_Probe-Matching.R** of Hutter et al., 2019)
- **03-2_Data-subsetting_JLW.R** (This is an extra step not in Hutter et al., 2019)
- **04_Loci_alignment_1May2019.R**
- **05_mtgenome_assembly_May8.R**
- **06_Trim_Align_Aug14.R**
- **07_Concat_CompleteMatrix_Aug29.R**
- **07-2_IQTREE_1May2019.R**
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











