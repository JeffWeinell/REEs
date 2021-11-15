#!/bin/bash
module load R

BLASTDBDIR=$1
VAR=$2
VAR2=$3

Rscript '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test3.R' "$BLASTDBDIR" "$VAR" "$VAR2"
