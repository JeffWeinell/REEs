# REEs Package

 - [Introduction](#Introduction)
 - [Dependencies](#Dependencies)
 - [Installation](#InstallingREEs)

<a name="Introduction"></a>
## Introduction.
This R package includes functions for identifying **rapidly-evolving exons (REEs)** using a set of whole genomes and at least one genome with exons annotated. I used some functions in this package to identify candidate phylogenomic target loci for targeted sequence capture in snakes ([SnakeCap](https://github.com/JeffWeinell/SnakeCap/blob/main/README.md)). These methods were inspired by [FrogCap](https://frogcap.com/) ([Hutter et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13517)) and RELEC ([Karin et al. 2019](https://doi.org/10.1093/molbev/msz263)).

<!--
The reasons for publishing these methods as an R package include (1) having a reproducible and citable pipeline for projects that use the SnakeCap probe set, and (2) to provide a method for researchers to select a set of loci for their phylogenomic studies.
-->

<a name="Dependencies"></a>
### Dependencies:
  - [R](https://www.r-project.org/) v3.6+
  - [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - [MAFFT](https://mafft.cbrc.jp/alignment/software/)

<a name="InstallingREEs"></a>
### Installation

To install REEs, use R package BiocManager and run the block of code below.

```
# Install and load BiocManager
install.packages(pkgs="BiocManager",repos = "http://cran.us.r-project.org")
library(BiocManager)

# Install REEs dependency packages
BiocManager::install(c("BSgenome","DECIPHER","phangorn","dplyr","data.table", "foreach","reutils","curl","knitr","devtools"),dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F)
BiocManager::install("gschofl/biofiles",dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F)

# Install REEs
BiocManager::install("JeffWeinell/REEs",update=FALSE,dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F,auth_token="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")
```

To install BLAST, follow the instructions [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/) or use the ```REEs::blast.install``` function.

```
# Load REEs
library(REEs)

# Install BLAST
blast.install()
```

To install MAFFT, follow the instructions [here](https://mafft.cbrc.jp/alignment/software/) and/or [here](https://mafft.cbrc.jp/alignment/software/installation_without_root.html) or use the ```REEs::install.mafft``` function.

```
# Load REEs
library(REEs)

# Install MAFFT
mafft.install()
```




