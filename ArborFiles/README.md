## Files provided to or by Arbor Biosciences:

**Version 1** targets and baits (not used; 120nt baits):
- Version1_Target-Loci_Jeff-Weinell_10Sep2018.fasta: First set of proposed target loci submitted to Arbor (corresponding to WeinellEntry1–WeinellEntry1898)
- WEINELL-ultrastringent-baits-180919je.fas.gz: Baits designed from version 1 targets
- WEINELL-ultrastringent-baitcoverage-180919je.txt.gz: bait coverage statistics for version 1 loci

**Version 2** targets and baits (not used; 120nt baits):
- Version1-loci-removed_ZeroBaitCoverageLoci.tsv: The list of version 1 targets removed after application of Arbor's ultrastringent filtering method. Removed loci included:
  - 70 REEs and 29 immune loci
- Version1-loci-removed_baits-nonspecific.tsv: list of targets filtered because all of their baits matched multiple regions of *T. sirtalis* genome. The loci removed for this reason included:
  - 8 REEs and 2 immune loci
- Version1-loci-removed_duplicate-targets.tsv: list of WeinellEntry targets filtered because they were already targeted under a different WeinellEntry name. Specifically, 12 immune loci were also identified as REEs and therefore these were initially included distinct targets.
- Version2_additional-targets_Entry1899to3152_20Sep2018.txt: New targets added to version 2, including:
  - 254 REEs: WeinellEntry1899–2152
  - 1,000 UCEs: WeinellEntry2153–3152
- Version 2 targets include 41 immune loci (plus 12 REEs that are also MHC loci), 1,990 REEs, and 1,000 UCEs.
  - ALLWEINELL-ultrastringent-baits-180925je.fas.gz: Baits designed from version 2 set of target loci
  - ALLWEINELL-ultrastringent-baits-180925je.fas.list.targcovg.table.gz: bait coverage statistics for version 2 loci.

**Version 3** targets and baits (not used; 120nt baits): 

- Version2-baits-removed_baits-nonspecific.txt: Probes from version 2 that were not included in version 3, because these hit multiple sites in the other genomes. **Need to figure out which genomes were blasted against in this step.**.
- Version2-loci-removed_baits-nonspecific.tsv: List of targets removed as a result of removing the multi-hit probes. The removed targets included:
  - 115 REEs, 14 immune loci, and 12 UCEs
- Version3_additional-targets_Entry3153to5735_4Oct2018.fasta: List of new targets added to version 3, including:
  - 2,025 short exon fragments: WeinellEntry3153–5177
  - 98 scalation loci: WeinellEntry5178-5275
  - 132 vision loci: WeinellEntry5276–5407
  - 328 ddRAD-like loci: WeinellEntry5408-5735.

- ADD2WEINELL-ultrastringent-baits-181005je.fas.gz: Baits designed from version 3 set of target loci

<!--
R Code used to get the list of loci in Version3-ZeroBaitCoverageLoci.tsv
library(ape)
loci.v3.added         <- read.dna(file="/Users/alyssaleinweber/Downloads/Weinell_Additional-Loci_Entry3153to5735_4Oct2018.fasta",format="fasta")
loci.v3.added.names   <- attributes(loci.v3.added)$names
probes.v3             <- read.dna(file="/Users/alyssaleinweber/Downloads/ADD2WEINELL-ultrastringent-baits-181005je.fas",format="fasta")
probes.v3.names       <- attributes(probes.v3)$dimnames[[1]]
probes.final          <- read.dna(file="/Users/alyssaleinweber/Downloads/Weinell_FinalProbeSet_20020Probes_7-Oct-2018.fasta",format="fasta")
probes.final.names    <- attributes(probes.final)$dimnames[[1]]
loci.v3.withProbes    <- unique(gsub("_.*","",probes.v3.names))
loci.final.withProbes <- unique(gsub("_.*","",probes.final.names))
loci.v3.zeroCoverage  <- setdiff(loci.v3.added.names,loci.v3.withProbes) ### these 74 loci were added to version 3, but had zero bait coverage, and therefore these were removed from version 4
-->

**Version 4** targets and baits (these are the baits we actually ordered and used for sequence capture; 120nt):

- Version3-loci-removed_ZeroBaitCoverageLoci.tsv: Version 3 targets removed by Arbor's ultrastringent filtering method; 74 targets were filtered, include 67 short exon fragments and 7 vision loci.
- Version3-baits-removed_others.tsv: List of version 3 baits not included in version 4 (n = 5,144), because the maximum number of baits is 20,020 for this kit.
- Version3-loci-removed_others.tsv: List of target loci removed (n = 2,189) as a result of removing the baits in Version3-probes-removed.tsv; removed targets included 212 REEs, 1,958 short exon fragments, 10 UCEs, 3 scalation, and 6 vision loci.
- Weinell_TargetLoci_Snakes_Final_18April2019.fa: This is the final set of target loci, for which baits were designed.
- Weinell_FinalProbeSet_20020Probes_7-Oct-2018.fasta: This is the set of baits purchased for sequence capture.
- Weinell_Final_SnakeCap_baits_target-coverage.txt: bait coverage table for the final set of target loci.

Description of Arbor's ultrastringent filtration:
"...  About our ultra-stringent filtration: a probe is eliminated from consideration if it is 25% RepeatMasked, or if its closest *Thamnophis* genomic hit is 25% or more soft-masked, or if the probe has multiple strong hybrid sites detected in the *Thamnophis* genome. To loosen the stringency we could increase the RM threshold and/or the number of tolerable hot hits in the genome. But I would strongly recommend you go with ultra-stringent if it hits sufficient target space for your application." - Arbor

**Future versions of target loci and baits** considering sequencing performance.

- Should use consensus sequences (across snakes) of target loci to design baits; this will reduce taxonomic dropout, which was not terrible but noticeable for version 4.
- Remove these loci, which were not captured for any species: 


