#!/bin/bash
# sbatch --nodes=1 --partition=sixhour --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/04.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupseqs/groups_fasta/" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs/" "NA"
module load anaconda
conda activate py36
INDIR=$1              # directory containing fasta files with unaligned sequences 
OUTDIR=$2             # directory where aligned sequences should be saved
GROUPSFILE=${3:-"NA"} # Optional file indicating which groups in $INDIR to evaluate
MODE=${4:-"A"}
SHij="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/04.sh"
PAR=250
MAXSEQS=100 # Alignment is skipped for input fasta files containing more sequences than MAXSEQS

# INDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupseqs/groups_fasta/"
# OUTDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupAlignments_SeqCapAndTargetRef_max100Seqs/"
# GROUPSFILE="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupseqs/TargetGroups_min4SeqCapSamples.txt" ### Not run with this

[ "$GROUPSFILE" = "NA" ] && INPATHS=$(find $INDIR -type f -name "group*") # ignores "ungrouped.fasta"
[ ! "$GROUPSFILE" = "NA" ] && INPATHS=$(awk '{ print "'$INDIR'"$1".fasta"}' $GROUPSFILE)
NUMPATHS=$(echo $INPATHS | wc -w)

[ ! -d "$OUTDIR" ] && mkdir -p "$OUTDIR"

if [ $MODE = "A" ] ; then
  for VARj in $(seq 1 $PAR $NUMPATHS) ; do
  #for VARj in $(seq 32501 $PAR 32750) ; do
    sbatch --nodes=1 --ntasks-per-node=10 --mem=10Gb --partition=sixhour $SHij $INDIR $OUTDIR $GROUPSFILE "B" $VARj
    # sbatch --nodes=1 --partition=sixhour $SHij $INDIR $OUTDIR $GROUPSFILE "B" $VARj
  done
fi

if [ $MODE = "B" ] ; then
  VARj=$5
  VARj2=$(($VARj+$PAR-1))
  [ $VARj2 -ge $NUMPATHS ] && VARj2=$NUMPATHS
  echo $(seq $VARj $PAR $VARj2) | awk '{print}'
  INPATHSi=$(echo $INPATHS | cut $"-f""$VARj""-""$VARj2" -d" ")
  NUMPATHSi=$(echo $INPATHSi | wc -w)
  for VARi in $(seq 1 $NUMPATHSi) ; do
    INPATHi=$(echo $INPATHSi | awk -v i=$VARi '{print $i}')
    NUMSEQSi=$(($(wc -l <$INPATHi) / 2 ))
    [ $NUMSEQSi -gt $MAXSEQS ] && continue
    OUTPATHi=$OUTDIR"/"$(basename $INPATHi)
    mafft --localpair --adjustdirection --maxiterate 1000 "$INPATHi" > "$OUTPATHi""_temp" &
  done
  wait
  for VARi in $(seq 1 $NUMPATHSi) ; do
    INPATHi=$(echo $INPATHSi | awk -v i=$VARi '{print $i}')
    OUTPATHi=$OUTDIR"/"$(basename $INPATHi)
    NUMSEQSi=$(( $(wc -l <$INPATHi) / 2 ))
    [ $NUMSEQSi -gt $MAXSEQS ] && continue
    seqkit seq -u -w 0 "$OUTPATHi""_temp" > "$OUTPATHi"
    rm "$OUTPATHi""_temp" &
  done
  wait
fi
