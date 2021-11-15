#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test18_iqtree.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes/" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/gene-trees/SeqCap_TrimmedGenomes_OneSequencePerSample/" "A"
module load anaconda
conda activate py36

### Captured variables
INDIR=$1    # Directory containing input alignments
OUTDIR=$2   # Output directory
MODE=${3:-"A"}
PAR=${4:-"10"}
SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test18_iqtree.sh"
###
# INDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes/" 
# OUTDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/gene-trees/SeqCap_TrimmedGenomes_OneSequencePerSample/"

[ ! -d "$OUTDIR" ] && mkdir -p "$OUTDIR"
LOGSDIR=$OUTDIR'/logfiles'
TREESDIR=$OUTDIR'/treefiles'
[ ! -d "$LOGSDIR" ] && mkdir -p "$LOGSDIR"
[ ! -d "$TREESDIR" ] && mkdir -p "$TREESDIR"

INPATHS=$(find $INDIR -type f | sort -V)
#### If some jobs fail, run the next two lines and then comment out the previous line and uncomment the line before NUMPATHS
## INPATHS=$(comm -23 <(ls "$INDIR" | sort) <(ls "$TREESDIR" | awk '{gsub(".treefile",".fasta")}1' | sort) | awk -v dir="$INDIR" '{print dir"/"$1}' | sort -V)
## echo "$INPATHS" > "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/filepaths_temp.txt"
## INPATHS=$(cat "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/filepaths_temp.txt")
NUMPATHS=$(echo "$INPATHS"| wc -l)

# SEQCAPNAMES=$WORKDIR"/samples_SeqCap.txt"

if [ $MODE = "A" ] ; then
  for VARj in $(seq 1 $PAR $NUMPATHS) ; do
    sbatch --nodes=1 --partition=sixhour $SHij $INDIR $OUTDIR "B" $VARj "$PAR"
  done
fi

### VARj slurm job
if [ $MODE = "B" ] ; then
  VARj=$4
  PAR=$5
  VARj2=$(($VARj+$PAR-1))
  [ $VARj -gt $NUMPATHS ] && exit
  [ $VARj2 -ge $NUMPATHS ] && VARj2=$NUMPATHS
  INPATHSi=$(echo $INPATHS | cut $"-f""$VARj""-""$VARj2" -d" ")
  OUTPATHSi=$(echo $INPATHSi | awk '{gsub(" ","\n")}1' | awk '{gsub("^.+[/]","")}1' | awk -v outdir=$OUTDIR"/" '{print outdir$1}')
  NUMPATHSi=$(echo $INPATHSi | wc -w)
  #echo "$INPATHSi" | awk '{gsub("^.*/","")}1'
  echo "$INPATHS" | awk -v X="$VARj" -v Y="$VARj2" 'NR>=X && NR<=Y' | awk '{gsub("^.*/","")}1'
  echo -e "\n\n"
  for VARi in $(seq 1 $NUMPATHSi) ; do
    INPATHi=$(echo $INPATHSi | awk -v i=$VARi '{print $i}')
    OUTPATHi=$(echo $OUTPATHSi | awk -v i=$VARi '{print $i}' | awk '{gsub(".fasta","")}1')
    TREEPATHi="$TREESDIR""/"$(echo $OUTPATHSi | awk -v i=$VARi '{print $i}' | awk '{gsub(".*/","")}1' | awk '{gsub(".fasta","")}1')
    [ -f "$TREEPATHi" ] && continue
    iqtree -s "$INPATHi" -pre "$OUTPATHi" -m TEST -bb 1000 &
  done
  wait
  ### Move output files to $LOGSDIR or $TREESDIR
  FILESPRESENT=$(find "$OUTDIR""/" -maxdepth 1 -type f)
  LOGFILES1=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".bionj"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  LOGFILES2=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".ckp.gz"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  LOGFILES3=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".contree"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  LOGFILES4=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".iqtree"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  LOGFILES5=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".log"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  LOGFILES6=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".mldist"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  LOGFILES7=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".model.gz"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  LOGFILES8=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".splits.nex"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  #LOGFILES9=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".uniqueseq.phy"}' | sort ) <(echo "$FILESPRESENT" | sort))
  TREEFILES=$(comm -12 <(echo "$OUTPATHSi" | awk '{gsub(".fasta","")}1' | awk '{print $1".treefile"}' | sort ) <(echo "$FILESPRESENT" | sort) )
  for LOGi in $(seq 1 "$(echo "$LOGFILES1" | wc -l)") ; do
    mv $(echo "$LOGFILES1" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  for LOGi in $(seq 1 "$(echo "$LOGFILES2" | wc -l)") ; do
    mv $(echo "$LOGFILES2" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  for LOGi in $(seq 1 "$(echo "$LOGFILES3" | wc -l)") ; do
    mv $(echo "$LOGFILES3" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  for LOGi in $(seq 1 "$(echo "$LOGFILES4" | wc -l)") ; do
    mv $(echo "$LOGFILES4" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  for LOGi in $(seq 1 "$(echo "$LOGFILES5" | wc -l)") ; do
    mv $(echo "$LOGFILES5" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  for LOGi in $(seq 1 "$(echo "$LOGFILES6" | wc -l)") ; do
    mv $(echo "$LOGFILES6" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  for LOGi in $(seq 1 "$(echo "$LOGFILES7" | wc -l)") ; do
    mv $(echo "$LOGFILES7" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  for LOGi in $(seq 1 "$(echo "$LOGFILES8" | wc -l)") ; do
    mv $(echo "$LOGFILES8" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  done
  wait
  #for LOGi in $(seq 1 "$(echo "$LOGFILES9" | wc -l)") ; do
  #  mv $(echo "$LOGFILES9" | awk -v i="$LOGi" 'NR==i') "$LOGSDIR" &
  #done
  #wait
  for TREEi in $(seq 1 "$(echo "$TREEFILES" | wc -l)") ; do
    mv $(echo "$TREEFILES" | awk -v i="$TREEi" 'NR==i') "$TREESDIR" &
  done
  wait
fi

