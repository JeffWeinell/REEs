
### Usage:

##### Load in SnakeCap functions and required R packages
```
source("~/SnakeCap_Functions.R") # load in the functions in this file
library(Biostrings)
library(stringr)
library(ape)
library(data.table)
```
##### Run make.partitioned.alignment separately for REEs, UCEs, ddRAD-like, MHC, scalation, and vision genes:

```
### Example usage for MHC genes:
make.partitioned.alignment(InputAlignmentFolder="~/Immune/unpartitioned/", output.dir="~/Immune/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")
```

