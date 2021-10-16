#!/bin/bash
# usage: sbatch --nodes=1 --ntasks-per-node=10 --mem=100Gb --time=6:00:00 --partition=sixhour '~/work/mitogenomes.sh' <"SampleName"> <"path/to/*READ1_unmerged.fastq.gz"> <"path/to/*READ2_unmerged.fastq.gz"> <"path/to/*_singletons.fastq.gz">
# NOTE: additional parameters must be set within "mitogenomes_script.R"

module load R
module load java
module load python

SNAME=$1
READ1=$2
READ2=${3:-"NULL"}
MERGED=${4:-"NULL"}

Rscript '/panfs/pfs.local/home/j926w878/work/mitogenomes_script.R' $SNAME $READ1 $READ2 $MERGED
