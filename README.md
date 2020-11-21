# SnakeCap Sequence Capture Probe Set:

## Contents

[Description of SnakeCap Probe Set](#Description)

[Methods](#Methods)
  - [Choosing loci to include in the probe set](#Methods.SelectingTargetLoci)
    - [Selecting the set of target Rapidly-evolving Exons (REEs)](#Methods.SelectingREEs)
    - [Selecting the set of target Ultraconserved-elements (UCEs)](#Methods.SelectingUCEs)
    - [Selecting the set of target ddRAD-like loci](#Methods.SelectingddRAD)
    - [Selecting the set of target functionally interesting (Functional) loci](#Methods.SelectingFunctional)
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
Functional: |  |  |  |  
Immune | (Exons of major histocompatibility complex (MHC) genes.) | 27 | 5,354 | 121–364
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

1. I downloaded the *Thamnophis sirtalis* genome, and its associated annotation table: **ref_Thamnophis_sirtalis-6.0_top_level_MHC.gff3** (n = 559,129 annotations). I renamed the sequences in the genome file to have the following format: **Thamnophis_sirtalis_GCF_001077635.1_read1**, etc., and wrote the renamed genome as a new file in sequential fasta format. I then made a two-column text file named **Scaffold-Name-Key.txt**, which has the following table (columns tab-delimited):

ScaffoldName | RefSeq.ScaffoldAccession
------------ | -------------
Thamnophis_sirtalis_GCF_001077635.1_read1 | RefSeq accession number (the ID from column 1 of the annotation table)
Thamnophis_sirtalis_GCF_001077635.1_read2 | ...
...|...
Thamnophis_sirtalis_GCF_001077635.1_readN|...

2. *Thamnophis sirtalis* exome extracted from genome. I used the function **filter.annotationTable** to: (1) filter the original annotation table to only include annotations for regions that are both CDS regions and ≥ 120bp in length (result: n = 115,907 annotations for 77,329 unique regions), and (2) write the filtered annotation table to a file: **CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_withGenBankAcc_longer120bp.gff3**. Then, I used the function **get.exome.from.annotationTable** to extract from the genome the DNA sequences included in the filtered annotation table, and to write the extracted DNA sequences (the exome) to a file in fasta format: **Thamnophis_sirtalis_exome_longer120bp.fas**.

3. Find *T. sirtalis* exons in other squamate genomes. I downloaded all squamate genomes available from NCBI (**Table 4**). I queried each *T. sirtalis* exon against each squamate genome using the program **tblastx** (allowing up to 50 matches to be saved per query). This step produces a Hit Table for each species, which includes stats on each query/target match, including the bitscore, which is a measurement how good the match is (how likely the match corresponds to homology). **Note to self**: I used the cluster submission file **TBLASTX.sh** to perform this step. Then, I filtered each of the full (i.e., 50 matches/query) hit tables to include only the best match/query (max bitscore) using the R function **reportBestMatches**.

4. Extract exomes (minimum exon length 120nt) for each squamate species. For each species and exon, I extracted the DNA sequence of the best match in the filtered hit table from step 2, and I saved these sequences to a fasta file; this step was performed using the function **get.exome.from.blastTable**.

5. Align shared exons and calculate stats. I used the R function **makeExomeStatsTable** to do all of the following:
  - obtain the set of *T. sirtalis* exons present in all exomes (from step 3), i.e., the shared exons
  - perform multiple sequence alignment (MAFFT algorithm) for each of the shared exons
  - calculate a set of stats for each shared exon alignment, and save results to the file **stats_exome_data_TBLASTX.txt** (results: includes stats for 66,489 alignments). 
  
Specifically, the following descriptors and stats were calculated for each exon alignment and saved to the file **stats_exome_data_TBLASTX.txt**: (1) NCBI Reference Sequence ID and location (stop_start) of CDS region on contig (2) number of species in alignment (always 11, because 11 species included, and only shared loci were aligned), (3) number of sites in which at least four species represented, (4) number of parsimony informative sites, (5) percent of sites parsimony informative, (6) mean pairwise percent genetic similarity to *T. sirtalis*, (7–17) percent genetic similarity to *T. sirtalis* for each species, (18) alignment width, (19) width of *T. sirtalis* exon, (20) gene name associated with exon, (21) mean number of variable sites compared to *T. sirtalis*, (22) minimum percent genetic similarity to *T. sirtalis* (among all species included), (23) minimum percent genetic similarity to *T. sirtalis* (among all snakes included). The gene name for each exon was extracted from the last column of the filtered annotation table.

6. Filter loci to include the most rapidly evolving loci that also meet criteria. I used the R function **pick.loci** to do the following: (1) filtered out loci if minimum percent genetic similarity (among snakes) to *T. sirtalis* was < 65% or = 100% (results = 64,546 loci retained; this temporary stats table was not written to a file). (2) For genes with multiple exons, I only kept the fastest evolving exon for each gene (the exon with the lowest mean pairwise genetic distance to *T. sirtalis*).

```
result.pident <- pick.loci(statsTable.path = "~/AlignedExonStats/stats_exome_data_TBLASTX.txt", output.dir = "~/AlignedExonStats/", primary.species = "Thamnophis_sirtalis", use.min.pident.subgroup = T, species.subgroup = c(7:14), min.pident.keep = c(65,100), max.capture.coverage = 1200000, write.stats.tables = T, plot.results = T, fast.stat = "pident")
```

Running the pick.loci function produced the table "result.pident", which contains percent identity stats for the 16,650 loci retained; this stats table was saved to the file  **stats_data_FastestExonPerGene.txt**.

From these 16,650 exons, I chose the set of exons that maximizes the number of variable sites while meeting the following constraints and conditions: total number of nucleotides targetted ≤ 1.2Mb and number of baits ≤ 20,000 (constraints of the 20K my-baits kit), bait size = 120nt, and bait tiling = 50% overlap.

```
Insert code showing the function used for selecting the optimal set of target exons given the constraints listed above. This script sorted the exons by % similarity to T. sirtalis (increasing simility), calculated the number of baits needed to capture each exon, and then performed a rolling summation across the vector of number of baits/exon, and a rolling summation of target exon length. The largest set of exons for which number of baits ≤ 20,000 and number of nucleotides ≤ 120,000,000 was chosen as the optimal set.
```

Result = 2,068 REEs retained; the stats table for these loci was written to the file **stats_data_FastestExonPerGene_best.txt**; an updated version of this stats table that includes the WeinellEntry locus names is **stats_data_FastestExonPerGene_best_20Nov2020.txt**.

These 2,068 REEs were submitted to Arbor Biosciences for probe design. Arbor performed ultrastringent filtration these loci which resulted in the removal of 70 REEs (**Version1_ZeroBaitCoverageLoci.tsv**), whereas 1,998 passed this step. An additional 123 REEs were removed because the baits designed to target these loci were all non-specific within the genomes of *T. sirtalis* and/or *Thermophis baileyi*; 212 other REEs were removed to allow some baits to be used to target other types of loci (UCEs, immune, scalation, vision, and ddRAD-like loci). The remaining 1,653 REEs were synthesized and included in the mybaits 20K bait kit (product no. 3001160). 


<!--
REEs that failed ultrastringent filtration: 70
REEs removed because all baits were non-specific: 123
REEs removed to make room for other loci: 212
Total REEs removed: 70+123+212 = 405
1,653 REEs were synthesized and included in the bait kit.
10 REEs unaccounted for...
-->

<!--
Version1_ZeroBaitCoverageLoci.tsv: 70 REEs, 29 immune loci
Version1_removed-loci_baits-nonspecific.tsv: 8 REEs, 2 immune loci
Version1_removed-loci_duplicate-targets.tsv: 12 immune loci
Version2-Loci-removed_allProbes-multiHit.tsv: 115 REEs, 14 immune loci, and 12 UCEs
Version3-loci-removed_ZeroBaitCoverageLoci.tsv: 67 short exon fragments, 7 vision loci
Version3-loci-removed_others.tsv: 212 REEs, 1,958 short exon fragments, 10 UCEs, 3 scalation, and 6 vision loci.
See the README file in ArborFiles folder for a description about how the bait kit changed after working with Arbor.
-->


<!--- Next step was done but wasn't necessary. It might be good to run this as a sanity check when picking REEs for future probe sets
7. I used the R function align.and.concatenate.best.exons to perform multiple sequence alignment (MAFFT algorithm) for each exon in the "stats_data_FastestExonPerGene_best.txt" file (which was generated by the pick.loci function in step 5), and to concatenate exon alignments and generate an associated partition file. Single-locus alignments, the concatenated loci alignment, and the partition file were each saved to file. I estimated gene trees from these alignments (using IQTREE) to get a sense for how much phylogenetic information each locus contained.
--->

<!--- Idk what I these commented out lines are for:
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

1. I downloaded the set of *Micurus fulvius* UCEs (n = 3,260) identified by Streicher and Wiens (2017). These were available as a fasta file called **micrurus_UCEs.fas**.

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

3. Then, I filtered the hit tables to include only the best match/query (max bitscore) using the R function **reportBestMatches**. Results were saved to the files: **Crotalus_horridus.blastn.UCEs.best.txt**, **Crotalus_mitchellii.blastn.UCEs.best.txt**, **Ophiophagus_hannah.blastn.UCEs.best.txt**, **Pantherophis_guttatus.blastn.UCEs.best.txt**, **Protobothrops_mucrosquamatus.blastn.UCEs.best.txt**, **Python_bivittatus.blastn.UCEs.best.txt**, **Thamnophis_sirtalis.blastn.UCEs.best.txt**, **Vipera_berus.blastn.UCEs.best.txt**.

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

<a name="Methods.SelectingFunctional"></a>
### Selecting the set of Functional loci

#### Overview:

Scalation loci: I targeted a subset of the genes included in the study by Holthaus et al. (2017). In that study, the authors identified homologous genes of the Epidermal Differentiation Complex (which are putatively involved in scalation) of *Python bivittatus* and *Ophiophagus hannah*. I downloaded the *Ophiophagus* scalation gene sequences using the table of genomic coordinates provided by Holthaus et al. (2017), and then used tblastn to search for and obtain homologous loci in *T. sirtalis*, *Protobothrops mucrosquamatus*, and *Crotalus horridus*.

Immune loci: I searched the annotation table of *T. sirtalis* for MHC (I or II) genes. Of these, only those for which Arbor Biosciences could create probes were included in the probe set.

Vision loci: I used blastn to search for the vision loci probes from Schott et al. (2017) (which were from *Anolis*, *Columba*, *Gallus*, and *Pelodiscus*, *Sceloporus*, or *Python*) within the snake genomes. Most of the SnakeCap probes for these loci are designed from *Ophiophagus* (n = 88), but some probes were designed from *Thamnophis* (n = 21), *Protobothrops* (n = 5), *Pantherophis* (n = 3), or *Python* (n = 2), when blastn of Schott et al 2017 probes did not yield a strong match in *Ophiophagus*.

#### Detailed, step-by-step methods for how I chose the set of functional loci:

1. 

2. 

3. 

<a name="ProbeSynthesis"></a>
### Probe Synthesis

After choosing the target loci, probes were designed by Arbor Biosciences with the following specifications: 50% tiling, 120nt/probe; 20,020 probes in total. See **Target-loci_Coverage_graph_22October2020.pdf** for a visual summary of target loci, probes, probe coverage, and features of loci including genes, mRNA/transcribed regions, and protein-coding (CDS) regions. This graph was generated with **graph_target_and_features.R** and then filesize reduction in Adobe Acrobat.

<a name="Sampling"></a>
### Taxa Sampled

<a name="LibraryPrep"></a>
### Sequence Capture Library Prep

Conducted by Arbor Biosciences; eight samples/pool; ...

<a name="DNASequencing"></a>
### DNA Sequencing

Novogene Illumina HiSeqX; paired-end sequencing with 150bp insert size.

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











