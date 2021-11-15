#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test10.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs/" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs_consensus/"
module load R
module load anaconda
conda activate py36

### RUN THIS CODE AFTER RUNNING MAFFT (04.sh) AND ONE OR MORE STEPS BEFORE RUNNING IQTREE ### 
# When mode="A", parallel subjobs will be created to obtain the consensus sequence for each input fasta in INDIR, and each consensus sequence will be written to a file in OUTDIR.
# If mode="C", OUTFILE will contain the consensus sequence for each input fasta in INDIR; consensus sequence names will reflect input fasta filenames.

### Captured variables
INDIR=$1  # Directory containing input alignments
OUT=$2   # Either the output directory (OUTDIR), modes A/B, or the output path for mode C.
MODE=${3:-"A"} # Mode A calls parallel subjobs, each of which have mode B. Mode C concatenates all output fasta files to OUTFILE, and should be run after the mode B jobs are all complete.
SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test10.sh"
Rij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test10.R"
###
## Modes A/B
# INDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs/"
# OUT="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs_consensus/"
## Modes 
# INDIRC="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs_consensus/"
# OUTC="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_consensus_InputSeqsMax100.fasta"

### Variables that sometimes need to be changed
INPATHS=$(find $INDIR -type f | sort -V)
NUMPATHS=$(ls $INDIR| wc -w)
PAR=500

if [ $MODE = "A" ] ; then
  OUTDIR="$OUT"
  [ ! -d "$OUTDIR" ] && mkdir -p "$OUTDIR"
  for VARj in $(seq 33501 $PAR $NUMPATHS) ; do
    sbatch --nodes=1 --ntasks-per-node=10 --mem=10Gb --partition=sixhour $SHij $INDIR $OUTDIR "B" $VARj
  done
fi

### VARj slurm job
if [ $MODE = "B" ] ; then
  OUTDIR="$OUT"
  VARj=$4
  VARj2=$(($VARj+$PAR-1))
  [ $VARj2 -ge $NUMPATHS ] && VARj2=$NUMPATHS
  # echo "$Rij"" ""$INDIR"" ""$OUTDIR"" ""$VARj"" ""$VARj2"
  Rscript "$Rij" "$INDIR" "$OUTDIR" "$VARj" "$VARj2"
fi

# After all VARj slurm jobs have completed
# 
# sbatch --nodes=1 --partition=sixhour "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test10.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs_consensus/" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_consensus_InputSeqsMax100.fasta" "C"
if [ $MODE = "C" ] ; then
  INDIRC=$1
  OUTC=$2
  cat $(find "$INDIRC" -type f) > $OUTC"_temp"
  seqkit seq $OUTC"_temp" -w 0 > "$OUTC"
  # OUTFILE="$VAR2"
  # seqkit scat -j 4 -I fasta -f $INDIR > "$OUTFILE"
fi



