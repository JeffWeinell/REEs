#!/bin/bash
# module load python/3.7
# module load perl/5.30
module load anaconda
conda activate py36

# usage # sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/RepeatMasker.sh" <"/path/to/input/sequence/file"> ["/path/to/output/directory"]

# Variables
SEQFILE="$1"
OUTDIR=${2:-PWD}

# Create output directory if it doesnt already exist
[ ! -d "$OUTDIR" ] && mkdir $OUTDIR

# Path to directory containing RepeatMasker Perl Script
RMDIR="/panfs/pfs.local/work/bi/bin/RepeatMasker/"
cd $RMDIR

# Running Perl Script
perl RepeatMasker -dir $OUTDIR $SEQFILE -species squamates

