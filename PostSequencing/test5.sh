#!/bin/bash
module load R

# Rscript '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5.R' "$1" "$2"
#VAR1=$1
#VAR2=$2
#VAR3=${3:-"1"}
#NUMFILES=$(ls -1 $INDIR | wc -l)
#VAR4=${4:-"1"}
Rscript '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5_v2.R' "$1" "$2" "$3" "$4"
