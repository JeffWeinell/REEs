# REEs Package

 - [Introduction](#Introduction)
 - [Software Dependencies](#Dependencies)
 - [Installing REEs](#InstallingREEs)
 - [Installing BLAST](#InstallingBLAST)
 - [Installing MAFFT](#InstallingMAFFT)
 - [How to use this package](#HowTo)
 - [Example](#Example)

<a name="Introduction"></a>
## Introduction.
This R package includes functions for identifying **rapidly-evolving exons (REEs)** using a set of whole genomes and at least one genome with exons annotated. I used some functions in this package to identify candidate phylogenomic target loci for targeted sequence capture in snakes ([SnakeCap])(https://github.com/JeffWeinell/SnakeCap/blob/main/README.md). These methods were inspired by [FrogCap](https://frogcap.com/) ([Hutter et al. 2021])(https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13517) and RELEC ([Karin et al. 2019])(https://doi.org/10.1093/molbev/msz263).

<!--
The reasons for publishing these methods as an R package include (1) having a reproducible and citable pipeline for projects that use the SnakeCap probe set, and (2) to provide a method for researchers to select a set of loci for their phylogenomic studies.
-->

<a name="Dependencies"></a>
### Software Dependencies:
  - [R](https://www.r-project.org/) v3.6+
  - [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - [MAFFT](https://mafft.cbrc.jp/alignment/software/)

<a name="InstallingREEs"></a>
### Installing REEs

Load R and then run the block of code below. Uncomment and modify the .libPaths() line to install REEs to a non-default location.
```
# Install and load BiocManager package
install.packages(pkgs="BiocManager",repos = "http://cran.us.r-project.org")
library(BiocManager)

# Install REEs dependencies
BiocManager::install(c("BSgenome","DECIPHER","phangorn","dplyr","data.table", "foreach","reutils","curl","knitr","devtools"),dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F)
BiocManager::install("gschofl/biofiles",dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F)

# Install REEs.
BiocManager::install("JeffWeinell/REEs",update=FALSE,dependencies=c("Depends", "Imports", "LinkingTo"),build_vignettes=F,auth_token="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")
```

<a name="InstallingBLAST"></a>
### Install BLAST

You can follow the instructions [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/) or use the REEs ```blast.install``` function to install BLAST.

```
# Load REEs
library(REEs)

# Install BLAST
blast.install()
```

<a name="InstallingMAFFT"></a>
### Install MAFFT

Instructions for installing MAFFT can be found [here](https://mafft.cbrc.jp/alignment/software/).
Also useful: [How to install MAFFT to a non-default location](https://mafft.cbrc.jp/alignment/software/installation_without_root.html).

Alternatively, use the ```install.mafft``` function.
```
# Load REEs
library(REEs)

# Install MAFFT
mafft.install()
```




