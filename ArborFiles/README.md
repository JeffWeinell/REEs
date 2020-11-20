## Files provided to or by Arbor Biosciences:

**Version 1** targets and baits (not used; 120nt baits):
- Target-Loci_Jeff-Weinell_10Sep2018.fasta: First set of proposed target loci submitted to Arbor (corresponding to WeinellEntry1–WeinellEntry1898)
- WEINELL-ultrastringent-baits-180919je.fas.gz: Baits designed from version 1 targets
- WEINELL-ultrastringent-baitcoverage-180919je.txt.gz: bait coverage statistics for version 1 loci

**Version 2** targets and baits (not used; 120nt baits):
- Version1_ZeroBaitCoverageLoci.tsv: The list of targets (version 1 set) filtered after application of Arbor's ultrastringent filtering method; these loci had 0% bait coverage.
- Other version 1 targets were filtered because they were identical to other positions in the genome, or because they were targeted twice (i.e., the positions were already targeted by another "WeinellEntry" target); this happened for some loci that were selected for inclusion for different reasons. Some of the immune and vision loci were already included because they were identified as single-copy rapidly-evolving exons (REEs). The following loci were filtered for one of these two afforementioned reasons: REEs: WeinellEntry44, WeinellEntry45, WeinellEntry496, WeinellEntry497, WeinellEntry589, WeinellEntry590, WeinellEntry1670, WeinellEntry1671; Immune loci: WeinellEntry1815, WeinellEntry1818, WeinellEntry1822, WeinellEntry1824, WeinellEntry1828, WeinellEntry1829, WeinellEntry1833, WeinellEntry1838, WeinellEntry1840, WeinellEntry1841, WeinellEntry1847, WeinellEntry1850, WeinellEntry1854, WeinellEntry1856, WeinellEntry1858, WeinellEntry1862, WeinellEntry1863, WeinellEntry1864, WeinellEntry1865, WeinellEntry1866, WeinellEntry1869, WeinellEntry1872, WeinellEntry1873, WeinellEntry1874, WeinellEntry1875, WeinellEntry1876. More details on these loci are shown in the following two tables.

Targets removed because sequence is non-specific in *Thamnophis* genome; only targets with non-zero bait coverage after ultrastringent filtering are shown here.
| loci having identical sequences  | *Thamnophis* scaffold and location | gene names used in genome annotation| type |
|---|---|---|---|
| WeinellEntry44; WeinellEntry45 | NW_013659646.1:15015-15975; NW_013659646.1:24066-25026 | vomeronasal type-2 receptor 26-like; vomeronasal type-2 receptor 116-like | REEs |
| WeinellEntry496; WeinellEntry497  | NW_013657804.1:817033-817753; NW_013657804.1:820482-821202 | olfactory receptor 10A2-like; olfactory receptor 10A2-like | REEs |
| WeinellEntry589; WeinellEntry590  | NW_013658610.1:34635-34875; NW_013658610.1:59791-60031 |  metalloproteinase inhibitor 4-like; metalloproteinase inhibitor 4-like | REEs |
| WeinellEntry1670; WeinellEntry1671  | NW_013658165.1:745869-746109; NW_013658165.1:748461-748701  | helicase with zinc finger domain 2-like; helicase with zinc finger domain 2-like  | REEs |
| WeinellEntry1815; WeinellEntry1820  | NW_013657739.1:619933-620293; NW_013657804.1:5056-5416 | class 2 MHC transactivator (exon7?); MHC class II transactivator-like | immune |

Target entries removed because already targeted by another entry:
| Duplicate pair target loci  | *Thamnophis* scaffold and location | locus filtered from version 2 targets | type (based on method/reason for including)|
|---|---|---|---|
| WeinellEntry397; WeinellEntry1818 | NW_013657739.1:625317-626997  | WeinellEntry1818 | REE; immune |
| WeinellEntry422; WeinellEntry1822 | NW_013657804.1:12345-14265   | WeinellEntry1822  | REE; immune |
| WeinellEntry559; WeinellEntry1824 | NW_013657921.1:388522-388882 | WeinellEntry1824  | REE; immune |
| WeinellEntry1539; WeinellEntry1828 | NW_013658177.1:704767-705247 | WeinellEntry1828 | REE; immune |
| WeinellEntry248; WeinellEntry1829 | NW_013658448.1:275425-275785 | WeinellEntry1829 | REE; immune |
| WeinellEntry612; WeinellEntry1833 | NW_013658549.1:435664-436024 | WeinellEntry1833  | REE; immune |
| WeinellEntry1232; WeinellEntry1838 | NW_013658556.1:485920-486280 | WeinellEntry1838  | REE; immune |
| WeinellEntry191; WeinellEntry1841 | NW_013659059.1:5414-5774 | WeinellEntry1841 | REE; immune |
| WeinellEntry836; WeinellEntry1847 | NW_013659625.1:112053-112413 | WeinellEntry1847 | REE; immune |
| WeinellEntry787; WeinellEntry1850 | NW_013659779.1:147829-148189 | WeinellEntry1850 | REE; immune |
| WeinellEntry728; WeinellEntry1854 | NW_013660884.1:27009-27369 | WeinellEntry1854 | REE; immune |
| WeinellEntry891; WeinellEntry1856 | NW_013661154.1:9614-9974 | WeinellEntry1856 | REE; immune |

- Weinell_Additional_Targets_20Sep2018.txt: Targets added to version 2 that were not in version 1, including more REEs (WeinellEntry1899–WeinellEntry2152) and 1,000 UCEs (WeinellEntry2153–WeinellEntry3152).
- Version 2 targets include 41 immune loci (plus 12 REEs that are also MHC loci), 1,990 REEs, and 1,000 UCEs.
- ALLWEINELL-ultrastringent-baits-180925je.fas.gz: Baits designed from version 2 set of target loci
- ALLWEINELL-ultrastringent-baits-180925je.fas.list.targcovg.table.gz: bait coverage statistics for version 2 loci

**Version 3** targets and baits (these are the baits we actually ordered and used for sequence capture; 120nt):

- Weinell_MultiHitProbes_Remove.txt: Probes from version 2 that were not included in version 3, because these hit multiple sites in the other genomes. **Need to figure out which genomes were blasted against in this step.**
- Version2-Loci-removed_allProbes-multiHit.tsv: List of targets removed as a result of removing the multi-hit probes. The removed targets included 115 REEs, 14 immune loci, and 12 UCEs.
- Weinell_Additional-Loci_Entry3153to5735_4Oct2018.fasta: Targets added to version 3 that were not in version 2; these included 2,025 **REEs?** (WeinellEntry3153–5177), 98 scalation loci (~WeinellEntry5178-WeinellEntry5275), 132 vision loci (WeinellEntry5276–5407), and 328 ddRAD-like loci (WeinellEntry5408-5735).
- ADD2WEINELL-ultrastringent-baits-181005je.fas.gz: Baits designed from version 3 set of target loci

**Version 4** targets and baits (these are the baits we actually ordered and used for sequence capture; 120nt):
- Version3-probes-removed.tsv: List of version 3 baits not included in version 4. (Use setdiff, gsub, and unique function in R to get this list by comparing the probe names in ADD2WEINELL-ultrastringent-baits-181005je.fas to the probe names in Weinell_FinalProbeSet_20020Probes_7-Oct-2018.fasta).
- Version3-Loci-removed.tsv: List of targets removed as a result of removing the baits in Version3-probes-removed.tsv.
- Weinell_FinalProbeSet_20020Probes_7-Oct-2018.fasta: This is the set of baits purchased for sequence capture.
- Arbor didn't provide a target coverage table for the final set of probes, but I calculated this myself and the results are in the file SnakeCap_probes_target-coverage.txt

"...  About our ultra-stringent filtration: a probe is eliminated from consideration if it is 25% RepeatMasked, or if its closest *Thamnophis* genomic hit is 25% or more soft-masked, or if the probe has multiple strong hybrid sites detected in the *Thamnophis* genome. To loosen the stringency we could increase the RM threshold and/or the number of tolerable hot hits in the genome. But I would strongly recommend you go with ultra-stringent if it hits sufficient target space for your application." - Arbor
