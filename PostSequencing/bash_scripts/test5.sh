#!/bin/bash
module load R
# Mode A/B/C usage: sbatch --nodes=1 --mem=100Gb --partition=sixhour --time=6:00:00 '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5.sh' "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupsMatrix_capVScap_capVStargets_v3.txt"
# Mode D usage: sbatch --nodes=1 --partition=sixhour "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5.sh"  "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupsMatrix_capVScap_capVStargets_v3.txt" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupseqs/N_UniqueSamples_Seqs_TargRefs_PerGroup_WithHeader_v2.txt" "D" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/allseqs_key2.txt"

# $1  # Input directory (mode A/B/C) or input path (mode D)
# $2  # Output path (all modes)
MODE=${3:-"A"}

# Number of group matrices to combine per subjob; only used for modes A/B
PAR=100
# Path to this sh file
SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5.sh"

### Mode A: Submits multiple parallel subjobs (each subjob has mode B) to concatenate the set of query i subject j group membership matrices into a new group membership matrix with groups recoded to avoid clashes between group indices.
# Run mode A/B/C before running mode D.
if [ $MODE = "A" ] ; then
  INDIR=$1
  OUTPATH=$2
  OUTDIRTEMP=$OUTPATH"_temp"
  NUMFILES=$(ls -1 $INDIR | wc -l)
  [ -d "$OUTDIRTEMP" ] && mkdir -p $OUTDIRTEMP
  for VAR in $(seq 1 $PAR $NUMFILES); do
    VARj=$VAR
    OUTPATHj=$OUTDIRTEMP"/groupsMatrix_"$VARj"to"$VARj2".txt"
    sbatch --nodes=1 --mem=100Gb --partition=sixhour --time=6:00:00 '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5.sh' "$INDIR" "$OUTPATHj" "B" "$VARj"
  done
  # wait
  # Need to add some code to check on all of the mode B subjobs.
  # sbatch --nodes=1 --mem=100Gb --partition=sixhour --time=6:00:00 '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5.sh' "NA" "$OUTPATH" "C" "$OUTDIRTEMP"
fi

if [ $MODE = "B" ] ; then
  INDIR=$1
  OUTPATHj=$2
  NUMFILES=$(ls -1 $INDIR | wc -l)
  VARj=$4
  VARj2=$(($VARj+$PAR-1))
  [ $VARj2 -ge $NUMFILES ] && VARj2=$NUMFILES
  Rscript '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5_v3.R' "$INDIR" "$OUTPATHj" "$VARj" "$VARj2" "$MODE" "NA" "NA"
fi

### Mode C: a subjob called by mode A after all mode B subjobs are submitted; not yet implemented.
if [ $MODE = "C" ] ; then
  OUTDIRTEMP=$4
  # While loop to periodically check on the status of mode B subjobs
  # If any mode B subjobs fail due to hitting a memory limit, submit new mode B subjobs for the associated VARj and VARj2.
  # When all mode B subjobs are complete, merge the partially merged matrices.
fi

### Mode D: generates a summary table for the the group membership matrix that was generated from mode A;
# Run mode D after running mode A/B/C.
# Each row of the summary table includes the following columns: (1) group,  (2) N_UniqueSamples, (3) N_Seqs, (4)  N_targetRef, (5) targetRef_seqNames, (6) targetRef_Names
# Note:  N_UniqueSamples includes the number of SeqCap plus targetRef samples. To get the number of unique SeqCap samples, subtract N_targetRef from N_UniqueSamples.
if [ $MODE = "D" ] ; then
  #### Need to pass the paths of the input matrices allseqsMat (output from...) and the groupsMat (output from mode A/B/C)
  GROUPSMATPATH=$1
  OUTPATH=$2
  ALLSEQSMATPATH=$4
  ### Must keep mode as the 5th variable of the R script
  Rscript '/panfs/pfs.local/home/j926w878/work/conda/snakecap/test5_v3.R' "$GROUPSMATPATH" "$OUTPATH" "NA" "$ALLSEQSMATPATH" "D" "NA" "NA"
fi

