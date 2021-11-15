#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test7_mergegroups.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupsMatrix_capVScap_capVStargets.txt"

### Captured variables
BLASTDBDIR=$1
GROUPMATPATH=$2
MODE=${3:-"A"}
SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test7_mergegroups.sh"
###
# BLASTDBDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10"
# GROUPMATPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupsMatrix_capVScap_capVStargets_v2.txt"
# GROUPMATPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupMatrix_consensusSeqCapAndTargetRef_vs_genomes.txt"

### Variables that sometimes need to be changed
GROUPSEQSDIR=$BLASTDBDIR"/groupseqs"
OUTDIR=$GROUPSEQSDIR"/groups_fasta"

## GROUPSEQSDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/"
## OUTDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/supergroups_fasta"

[ ! -d "$OUTDIR" ] && mkdir $OUTDIR

PAR=100
## MINSAMPLES=4
if [ $MODE = "A" ] ; then
  # Groups names (1, 2, ..., N; ungrouped)
  #GROUPNAMES=$(awk -F '\t' 'NR>1{print "group"$2}' $GROUPMATPATH | sort | uniq | sort -V)" ungrouped"
  GROUPNAMES=$(awk -F '\t' 'NR>1{print $2}' $GROUPMATPATH | sort | uniq | sort -V)
  # For each group, gather all sequences of the group and put them in a fasta file.
  NGROUPS=$(echo $GROUPNAMES | wc -w)
  LASTGROUP=$(echo $GROUPNAMES | awk -v i=$NGROUPS '{print $i}')
  #GROUPSETS=$(seq 1 $PAR $NGROUPS)
  #NSETS=$(echo $GROUPSETS | wc -w)
  #LASTSETMAX=$(($PAR*$NSETS))
  touch $GROUPSEQSDIR"/groupsets.txt"
  COUNTER=0
  for VARj in $(seq 1 $PAR $NGROUPS) ; do
    COUNTER=$(($COUNTER+1))
    # GROUPNAMEi=$(echo $GROUPNAMES | awk -v i=$VARj '{print $i}')
    # declare -i VARj2=$VARj+$PAR-1
    VARj2=$(($VARj+$PAR-1))
    [ $VARj2 -ge $NGROUPS ] && VARj2=$NGROUPS
    # $(seq $VARj1 $VARj2)
    # SETi=$(echo $GROUPNAMES | cut -f1-100 -d" ")
    SETj=$(echo $GROUPNAMES | cut $"-f""$VARj""-""$VARj2" -d" ")
    echo $SETj >> $GROUPSEQSDIR"/groupsets.txt"
    # echo $(seq 1 50) | awk '{print "group"$0".fasta"}' | awk '{gsub(" ",".fasta group")}1'
    # cat $(find $GROUPSEQSDIR -regex '.*[^/]*_'${GROUPNAMEi}'.fasta') > $GROUPSEQSDIR"/"$GROUPNAMEi".fasta"
    ## cat $(find $GROUPSEQSDIR -regex '.*[^/]*_'${GROUPNAMEi}'.fasta') > $OUTDIR"/"$GROUPNAMEi".fasta"
    sbatch --nodes=1 --partition=sixhour $SHij $BLASTDBDIR $GROUPMATPATH "B" $COUNTER
    #for VARi in $(seq 1 $PAR) ; do
    #  GROUPNAMEij=$(echo $SETj | awk -v i=$VARi '{print $i}')
    #  cat $(find $GROUPSEQSDIR -regex '.*[^/]*_'${GROUPNAMEij}'.fasta') > $OUTDIR"/"$GROUPNAMEij".fasta" &
    #done
  done
fi

if [ $MODE = "B" ] ; then
  COUNTER=$4
  SETj=$(sed "${COUNTER}q;d" $GROUPSEQSDIR"/groupsets.txt")
  for VARi in $(seq 1 $PAR) ; do
    GROUPNAMEij=$(echo $SETj | awk -v i=$VARi '{print $i}')
    cat $(find $GROUPSEQSDIR -regex '.*[^/]*_'${GROUPNAMEij}'.fasta') > $OUTDIR"/"$GROUPNAMEij".fasta" &
  done
  wait
fi






