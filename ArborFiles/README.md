## Files provided to or by Arbor Biosciences:

**Version 1** targets and baits (not the final set; 120nt baits):
- Target-Loci_Jeff-Weinell_10Sep2018.fasta: First set of proposed target loci submitted to Arbor (corresponding to WeinellEntry1–WeinellEntry1898)
- WEINELL-ultrastringent-baits-180919je.fas.gz: Baits designed from version 1 targets
- WEINELL-ultrastringent-baitcoverage-180919je.txt.gz: bait coverage statistics for version 1 loci

**Version 2** targets and baits (not the final set; 120nt baits):
- Targets and baits in version 1 that were filtered/not included in version 2: all loci with 0% bait coverage from ultrastringent filtering of version 1 , plus the following loci: WeinellEntry44, WeinellEntry45, WeinellEntry496, WeinellEntry497, WeinellEntry589, WeinellEntry590, WeinellEntry1670, WeinellEntry1671, WeinellEntry1815, WeinellEntry1818, WeinellEntry1822, WeinellEntry1824, WeinellEntry1828, WeinellEntry1829, WeinellEntry1833, WeinellEntry1838, WeinellEntry1840, WeinellEntry1841, WeinellEntry1847, WeinellEntry1850, WeinellEntry1854, WeinellEntry1856, WeinellEntry1858, WeinellEntry1862, WeinellEntry1863, WeinellEntry1864, WeinellEntry1865, WeinellEntry1866, WeinellEntry1869, WeinellEntry1872, WeinellEntry1873, WeinellEntry1874, WeinellEntry1875, WeinellEntry1876.
- Targets and baits added to version 2 that were not in version 1: 
- ALLWEINELL-ultrastringent-baits-180925je.fas.gz: WeinellEntry1899–WeinellEntry2152 and WeinellEntry2153–WeinellEntry3152.
- ALLWEINELL-ultrastringent-baits-180925je.fas.list.targcovg.table.gz: bait coverage statistics for version 2 loci

**Version 3** targets and baits (these are the baits we actually ordered and used for sequence capture; 120nt):
- 
- ADD2WEINELL-ultrastringent-baits-181005je.fas.gz: This is the set of baits purchased for sequence capture.
- Arbor didn't provide a target coverage table for the final set of probes, but I calculated this myself and the results are in the file SnakeCap_probes_target-coverage.txt

"...  About our ultra-stringent filtration: a probe is eliminated from consideration if it is 25% RepeatMasked, or if its closest *Thamnophis* genomic hit is 25% or more soft-masked, or if the probe has multiple strong hybrid sites detected in the *Thamnophis* genome. To loosen the stringency we could increase the RM threshold and/or the number of tolerable hot hits in the genome. But I would strongly recommend you go with ultra-stringent if it hits sufficient target space for your application." - Arbor
