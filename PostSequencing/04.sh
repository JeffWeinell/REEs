#!/bin/bash
module load R/4.0
TARGETS=$1
WORKDIR=${2:-$PWD}
ISTART=${3:-$(echo 1)}
IEND=$4
Rscript '/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/04_Loci_alignment_23Sep2022_SnakeCap2.R' $TARGETS $WORKDIR $ISTART $IEND
