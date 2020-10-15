# SnakeCap Sequence Capture Probe Set:

### Description of SnakeCap Probe Set

The probe set includes 3,129 single-copy loci (1,517,011 nt) shared across snakes. Loci are categorized into four types: (1) rapidly evolving exons (REEs; n = 1,653), (2) ultra-conserved elements (UCEs; n = 907), (3) ddRAD-like loci (n = 328), (4) and functionally interesting genes, which includes 27 major histocompatibility complex (MHC) genes, 119 vision genes, and 95 scalation genes.

REEs include one or more entire exons and one or both exon-flanking regions, and range in length from 121 to 7,501 nt. I used a modified version of the FrogCap pipeline (Hutter et al., 2019) to select the optimal set of REEs from an alignment of snake exomes.

SnakeCap UCEs are a subset of the Micrurus fulvius UCEs from Streicher et al. (2017).

ddRAD-like loci are shared, single-copy loci identified from in-silico ddRAD using recognition sites for SbfI and EcoRI restriction enzymes.

Functional loci included entire or partial gene regions that have previously been predicted or known to function in either (1) vertebrate immune systems, (2) vision, (3) or scalation.

### Methods

#### Selecting the set of target REEs

##### Overview: 

##### Detailed, step-by-step methods for how I chose the set of target REEs

1. I downloaded the Thamnophis sirtalis genome and the associated annotation table ("ref_Thamnophis_sirtalis-6.0_top_level_MHC.gff3" [n = 559,129 annotations]). I renamed the sequences in the T. sirtalis genome (e.g.: "Thamnophis_sirtalis_GCF_001077635.1_read1") and wrote the genome as a new file in sequential fasta format. I then made a two-column text file named: "Scaffold-Name-Key.txt", which has the new genome read names in the first column (column name = "ScaffoldName"), and the RefSeq accession number (i.e., the ID from column 1 of annotation table) the in the second column (column name = "RefSeq.ScaffoldAccession"). I may have used a code from Carl Hutter to generate the sequential fasta genome with renamed sequences.

2. Thamnophis sirtalis exome extracted from genome. I used the function I used the function filter.annotationTable to: (1) filter the original annotation table to only include annotations for regions that are both CDS regions and ≥ 120bp in length (n = 115,907 annotations; 77,329 unique regions), and (2) write the filtered annotation table to a file ("CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_withGenBankAcc_longer120bp.gff3"). Then, I used the function get.exome.from.annotationTable to extract from the genome the DNA sequences included in the filtered annotation table, and to write the extracted DNA sequences (i.e. exome) to a file in fasta format ("Thamnophis_sirtalis_exome_longer120bp.fas ").

3. Find T. sirtalis exons in other squamate genomes. I downloaded all squamate genomes available from NCBI (Table 4). I queried each T. sirtalis exon (from step 1) against each squamate genome using the program tblastx (up to 50 matches saved per query). This step produces a Hit Table for each species that includes stats on each query/target match, including the bitscore, which is a measurement how good the match is (i.e., how likely the match corresponds to homology). I used the cluster submission file TBLASTX.sh to perform this step. Then, I filtered each of the full (i.e., 50 matches/query) hit tables to include only the best match/query (max bitscore) using the R function reportBestMatches.

4. Extract exomes (minimum exon length 120nt) for each squamate species. For each species, I extracted the DNA sequences corresponding to each best match in the filtered hit table (from step 2), and I saved these sequences (i.e. exome) to a fasta file. This step was performed using the function get.exome.from.blastTable.

5. Align shared exons and calculate stats. I used the R function makeExomeStatsTable to do all of the following: (1) obtain the set of T. sirtalis exons present in all exomes (from step 3), i.e., the shared exons; (2) perform multiple sequence alignment (MAFFT algorithm) for each of the shared exons; (3) calculate a set of stats for each shared exon alignment, and save results to the file "stats_exome_data_TBLASTX.txt" (results: includes stats for 66,489 alignments). Specifically, the following descriptors and stats were calculated for each exon alignment: (1) NCBI Reference Sequence ID and location (stop_start) of CDS region on contig (2) number of species in alignment (always 11, because 11 species included, and only shared loci were aligned), (3) number of sites in which at least four species represented, (4) number of parsimony informative sites, (5) percent of sites parsimony informative, (6) mean pairwise percent genetic similarity to T. sirtalis, (7–17) percent genetic similarity to T. sirtalis for each species, (18) alignment width, (19) width of T. sirtalis exon, (20) gene name associated with exon, (21) mean number of variable sites compared to T. sirtalis, (22) minimum percent genetic similarity to T. sirtalis (among all species included), (23) minimum percent genetic similarity to T. sirtalis (among all snakes included). The gene name for each exon was extracted from the last column of the filtered annotation table. The stats table was saved to a file named "stats_exome_data_TBLASTX.txt".

6. Filter loci to include the most rapidly evolving loci that also meet criteria. I used the R function pick.loci to do the following: (1) filtered out loci if minimum percent genetic similarity (among snakes) to Thamnophis sirtalis was < 65% or =100% (results = 64,546 loci retained; this stats table was not written to a file). For genes with multiple exons, I only kept the fastest evolving exon for each gene (i.e., only the exon with the lowest mean pairwise genetic distance to T. sirtalis; results = 16,650 loci retained; stats table: "stats_data_FastestExonPerGene.txt"). Next, I retained the set of exons with maximum number of variable sites in the target loci set, given a constraint on the total number of nucleotides that can be targeted (=1.2Mb for the SnakeCap probe set designed from the 20K Mybaits Kit; results = 2,068 loci retained [actually 2,071 because initially used max.capture.coverage = 1201000]; stats table was written to the file "stats_data_FastestExonPerGene_best.txt").

```

result.pident <- pick.loci(statsTable.path = "~/AlignedExonStats/stats_exome_data_TBLASTX.txt", output.dir = "~/AlignedExonStats/", primary.species = "Thamnophis_sirtalis", use.min.pident.subgroup = T, species.subgroup = c(7:14), min.pident.keep = c(65,100), max.capture.coverage = 1200000, write.stats.tables = T, plot.results = T, fast.stat = "pident")

```

7. (optional). I used the R function align.and.concatenate.best.exons to perform multiple sequence alignment (MAFFT algorithm) for each exon in the "stats_data_FastestExonPerGene_best.txt" file (which was generated by the pick.loci function in step 5), and to concatenate exon alignments and generate an associated partition file. Single-locus alignments, the concatenated loci alignment, and the partition file were each saved to file. I estimated gene trees from these alignments (using IQTREE) to get a sense for how much phylogenetic information each locus contained.

8. For each locus, choose target region that will be used for probe design, considering the length of each target exon, probe length (120bp), and tiling regime (50% overlap between adjacent probes).

Li = length of ith target exon
nprobes,i = minimum # of probes needed to cover Li
Lprobes,i = # of nucleotides covered by nprobes,i
nflank,i = # of exon-flanking nucleotides targeted (total, upstream + downstream)
nflank.5',i = # of nucleotides targeted upstream of exon
= (nflank,i)/2 rounded down to nearest integer
nflank.3',i = # of nucleotides targeted downstream of exon
= nflank,i - nflank.5',i


#### Selecting the set of target UCEs

##### Overview:

Target UCEs include 907 of the 3,260 UCEs previously identified in Micrurus fulvius (Streicher et al., 2017; Table X). First, I filtered the full set of Micrus UCEs to only include those present in all NCBI snake genomes (n = 2,968 UCEs). From this shared set of UCEs I (probably; check script to confirm) selected the 1,000 loci with the greatest number of informative sites, or, that were on different T. sirtalis contigs... Of these, 907 were included in the final probe set (following Arbor Biosciences' protocols for probe filtering).

##### Detailed, step-by-step methods for how I chose the set of target UCEs:

1. I downloaded the set of Micurus fulvius UCEs (n = 3,260) identified by Streicher et al. (2017). These were available as a fasta file, with the name "micrurus_UCEs.fas".

2. I searched for Micrurus UCEs in the other snake genomes. I queried each Micrurus UCE against each squamate genome using the NCBI blastn algorithm (saving ≤ 50 matches per query), which was implemented using the R wrapper function blastnR. Then, I filtered each of the 50-match hit tables to include only the best match/query (max bitscore) using the R function reportBestMatches.

3. From each genome, extract the best-match UCEs. From each snake genome, I extracted the DNA sequences corresponding to each best match in the filtered UCE hit table, and I saved these sequences to a fasta file. This step was performed using the function get.UCEs.from.blastTable.

4. Align the set of UCEs found in all snake genomes. This was done using the R function align.bestHit.UCEs.

5. I selected 1,000 UCEs subset of UCEs that... were on different T. sirtalis contigs, or, that had the most phylogenetic information or the largest mean pairwise genetic distance. I need to check on this...

#### Selecting the set of target ddRAD-like loci

##### Overview:

##### Detailed, step-by-step methods for how I chose the set of target ddRAD-like loci:

1. Search (blast?) for Sbfi recognition site in T. baileyi genome (sense strand contigs); output = a NCBI-format hit table
2. Search (blast?) for EcoRI recognition site in T. baileyi genome (sense strand contigs); output = a NCBI-format hit table
3. Search (blast?) for Sbfi recognition site in T. baileyi genome (antisense strand contigs); output = a NCBI-format hit table
4. Search (blast?) for EcoRI recognition site in T. baileyi genome (antisense strand contigs); output = a NCBI-format hit table
5. Filtered T. baileyi contigs (sense strand) to only include those with both restriction enzyme recognition sites.
6. Filtered T. baileyi contigs (antisense strand) to only include those with both restriction enzyme recognition sites.
7. For the set of contigs containing both recognition sites, extract the region between each pairwise combination of RE sites.
8. Filter extracted regions to keep only those with length between 900–1000bp.
9. BLAST (tblastx, tblastn, blastx, blastn?) each sequence in the set of 900-1000bp extracted regions to search within each snake genome
10. Keep the set of single-copy sequences present in all snakes genomes, and design probes for these target loci.

Calculated the pairwise genomic distance distance between 
Set of 900–1000bp regions of the Sense Strand containing Sbfi and EcoRI recognition sites: "ddRAD-like-loci_SenseStrand_SbfI-EcoRI_900to1000bp_PASSED_HitTable.txt"

#### Selecting the set of Functional loci

##### Overview:

Scalation loci: I targeted a subset of the genes included in the study by Holthaus et al. (2017). In that study, the authors identified homologous genes of the Epidermal Differentiation Complex (which are putatively involved in scalation) of Python bivittatus and Ophiophagus hannah. I downloaded the Ophiophagus scalation gene sequences using the table of genomic coordinates provided by Holthaus et al. (2017), and then used tblastn to search for and obtain homologous loci in Thamnophis sirtalis, Protobothrops mucrosquamatus, and Crotalus horridus.

Immune loci: I searched the annotation table of Thamnophis sirtalis for MHC (I or II) genes. Of these, only those for which Arbor Biosciences could create probes were included in the probe set.

Vision loci: 

##### Detailed, step-by-step methods for how I chose the set of functional loci:

1.

2.

3.



## References:

Hutter C.R., Cobb K.A., Portik D., Travers S., Wood Jr. P.L., and Brown R.M. (2019). FrogCap: A modular sequence capture probe set for phylogenomics and population genetics for Anurans, assessed across multiple phylogenetic scales. bioRxiv 825307. doi: https://doi.org/10.1101/825307.

Streicher et al. (2017).

Holthaus et al. (2017).


