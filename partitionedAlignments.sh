#!/bin/bash
# usage: sbatch --nodes=1 --ntasks-per-node=30 --time=6:00:00 --partition=sixhour "/panfs/pfs.local/home/j926w878/work/conda/snakecap/partitionedAlignments.sh"
module load R

### Description:
# Runs the REEs function make.partitioned.alignment()

STARTi=$1
ENDi=$2

BLASTDBDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/SnakeCap_AllSamples/blast_database11/"
INPATH=$BLASTDBDIR"/alignments_WENames"
OUTPATH=$BLASTDBDIR"/partitionedAlignments"
TARGETS_CDS_PATH="/panfs/pfs.local/home/j926w878/work/SequenceCapture/Weinell_TargetLoci_Snakes_Final_targetCDS_v4.fa"
SAMPLESKEYPATH=$BLASTDBDIR"/allsamples_key.txt"
NCORES=$(nproc --all)
Rscript <(echo -e "
.libPaths('/panfs/pfs.local/home/j926w878/work/R-packages')
library(REEs)
INPATH <- '"$INPATH"'
OUTDIR <- '"$OUTPATH"'
TARGETS_CDS_PATH <- '"$TARGETS_CDS_PATH"'
allsampleskey <- read.csv('"$SAMPLESKEYPATH"',sep="\t",header=F)
X1="$STARTi"
X2="$ENDi"
THREADS=max(floor(($NCORES)/2),1)
locusi  <- REEs::make.partitioned.alignment(input.path=INPATH,output.dir=OUTDIR,TargetCDS.path=TARGETS_CDS_PATH,ith.locus.start=X1,ith.locus.end=X2,mafft.params=sprintf('--auto --adjustdirection --op 3 --ep 0.123 --quiet --thread %s',THREADS),trimto=allsampleskey[,'V2'][allsampleskey['V3']=='SeqCap'])
")


