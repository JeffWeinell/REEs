#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --mem=50Gb --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test15_cropToSequence.sh"
# module load R
module load anaconda
conda activate py36

### Description:
# 
# 
# Note: dont unalign seqs in this script because still need to flatten sliding windows sequences into a single sequence; need something like paired-end read merging but for sliding windows
# 

# Path to this bash script. 
# SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test15_cropToSequence.sh"

# UADIR=$1  # Path to directory containing input, untrimmed alignments
# TADIR=$2  # Path to directory where output, trimmed (cropped) alignments should be saved; sequences include ConsensusSCTR and genome samples
# TGDIR=$3  # Path to directory where output, trimmed (cropped) unaligned (gaps-removed) sequences of genomic samples should be saved

#MINSEQLENGTH=100
#PAR=1

### untrimmed alignments
UADIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCapTargRefConsensus_GenomesUntrimmed_max100seqs/"
#UADIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCapTargRefConsensus_GenomesUntrimmed_max100seqs_MAFFT-auto-windows"
### trimmed alignments
TADIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCapTargRefConsensus_GenomeSeqs_TrimmedToConsensus"
# TGDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroups_GenomeSeqs_TrimmedToConsensus"
STATSMATPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_ConsensusSCTR_vs_GenomesUntrimmed_STATS.txt"
#SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test15_cropToSequence.sh"
# path to supergroups defined in test12.sh
GroupStringsPath_OLD="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/GroupsStrings_Consensus-vs-Genomes_ConsensusSeqs.txt"
# path to the output file that will contain the subset of lines from $GroupStringsPath_OLD that define supergroups that were successfully aligned/trimmed; the other lines hereafter do not define supergroups
GroupStringsPath_NEW="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupStrings_groups.txt"
#STATSMATPATH_TA="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_ConsensusSCTR_vs_GenomesTrimmed_STATS.txt"
#STATSMAT_TA=$(seqkit stats $(find "$TADIR" -type f -name *fasta | sort -V ) -e -b )
#echo "$STATSMAT_TA" > "$STATSMATPATH_TA"

[ ! -d "$TADIR" ] && mkdir $TADIR
#[ ! -d "$TGDIR" ] && mkdir $TGDIR

# V1=${1:-"NA"}
# V2=${2:-"NA"}

V1=$1
V2=$2

## if [ "$V1"="NA" ] && [ "$V2"="NA" ] ;
##    then
##    # Create the STATS file summarizing the alignments that were generated from running test14.sh
##    STATSMAT=$(seqkit stats $(find "$UADIR" -type f -name *fasta | sort -V ) -e -b )
##    echo "$STATSMAT" > "$STATSMATPATH"
##    # Exit with zero exit status
##    exit
## fi

#echo "test"
#exit

## if [ -f "$STATSMATPATH" ] ; 
##   then 
##   STATSMAT=$(cat "$STATSMATPATH") ; 
##   else 
##   STATSMAT=$(seqkit stats $(find "$UADIR" -type f -name *fasta | sort -V ) -e -b )
##   echo "$STATSMAT" > "$STATSMATPATH"
## fi
## 

STATSMAT=$(cat "$STATSMATPATH")

# FILES=$(seqkit stats $(find "$UADIR" -type f -name *fasta | sort -V ) -e -b | awk 'NR>1{if($5 > 0) {print $1}}')
#FILES_OG=$( echo "$STATSMAT" | awk 'NR>1{if($5 > 0) {print $1}}')
#FILES=$(comm -3 <(echo "$FILES_OG" | sort) <(ls "$TADIR" | sort) | sort -V)
FILES=$( echo "$STATSMAT" | awk 'NR>1{if($5 > 0) {print $1}}')
NUMFILES=$(echo "$FILES" | wc -l)

[ "$V2" -gt "$NUMFILES" ] && V2="$NUMFILES"

for VARi in $(seq "$V1" "$V2"); do
   ## Filepath to current untrimmed supegroup alignment
   #ALNUPATHi=$(echo "$ALNUPATHS" | awk -v i="$VARi" '{print $i}')
   ## Basename of $ALNUPATH
   #FILENAME=$(echo "$ALNUPATHi" | awk '{gsub("^.+/","")}1')
   FILEi=$(echo $FILES | awk -v i=$VARi '{print $i}')
   ALNUPATHi="$UADIR""/""$FILEi"
   # Read in the untrimmed alignment
   ALNU=$(cat "$ALNUPATHi")
   # Get the ConsensusSCTR sequence(s) from $ALNU
   ALNU_CON=$( echo "$ALNU" | seqkit grep --use-regexp --pattern 'consensus')
   # get the alignment position of the first and last non-gap character of each ConsensusSCTR sequence
   RANGES_CON=$(echo "$ALNU_CON" | seqkit locate --ignore-case --pattern "[A-Z].+[A-Z]" --use-regexp --hide-matched | awk '{if($4=="+"){print $5,$6}}')
   # get the alignment position of the first and last non-gap character among all of the ConsensusSCTR sequences in the supergroup alignment
   FIRST_CON=$(echo "$RANGES_CON" | awk '{print $1"\n"$2}' | sort -V | head -n 1)
   LAST_CON=$(echo "$RANGES_CON" | awk '{print $1"\n"$2}' | sort -V | tail -n 1)
   # Extract the region from $FIRST_CON to $LAST_CON of $ALNU, and then drop sequences with only gaps
   ALNT=$(echo "$ALNU" | seqkit subseq --region "$FIRST_CON":"$LAST_CON" | seqkit grep --by-seq --invert-match --use-regexp --pattern ^-*$)
   # Get genome sample sequences from $ALNT and remove gaps to unalign
   # SEQS_GEN=$(echo "$ALNT" | seqkit grep -r -p 'sample' | seqkit replace -p "-" -s)
#   SEQS_GEN=$(echo "$ALNT" | seqkit grep --use-regexp --pattern 'sample' | seqkit seq --remove-gaps --min-len 1)
   # Skip to next alignment if all sequences are less than $MINSEQLENGTH nongap characters.
#   SEQS_PASS=$(echo "$SEQS_GEN" | seqkit fx2tab --length --name | awk -v minlen="$MINSEQLENGTH" '{if($2>=minlen){print $1}}')
#   N_SEQS_PASS=$(echo "$SEQS_PASS" | wc -l)
   N_SEQS_PASS=$(echo "$ALNT" | seqkit stats | awk 'NR==2{print $4}')
   [ "$N_SEQS_PASS" -eq 0 ] && continue
   # Save $ALNT to $TADIR, and use the same filename as that of $ALNUPATH
   echo "$ALNT" | seqkit seq -w 0 > "$TADIR""/""$FILEi"
   # Save $SEQS_GEN to $TGDIR, and use the same filename as that of $ALNUPATH
#   echo "$SEQS_GEN" | seqkit seq -w 0 > "$TGDIR""/""$FILEi"
#   STATS1=$(seqkit stats "$TADIR""/""$FILEi" -e -b)
#   STATS2=$(seqkit stats "$TGDIR""/""$FILEi" -e -b | awk 'NR==2')
#   STATS=$(cat <(echo "$STATS1") <(echo "$STATS2"))
#   echo "$STATS"
done

FILESi=$(echo "$STATSMAT" | awk 'NR>1 {if($5 > 0){print $1}}' | awk -v x="$V1" -v y="$V2" 'NR>=x && NR<=y{print $1}')
OUTPATHSi=$(echo "$FILESi" | awk -v dir="$TADIR" '{print dir"/"$1}')
seqkit stats $OUTPATHSi -b -e

###### Define supergroup strings for the supergroups that could be aligned; those that could not be aligned are not longer valid, and any groups that were contained in those unalignable "supergroups" can be considered as their own "supergroup" (non-genome containing) beginning at step 16.
## Supergroups sensu $GroupStringsPath 
OLD_SG=$(awk -F'|' 'NF=1' "$GroupStringsPath_OLD")
## Alignments "$TADIR"
PASSED_SG=$(find "$TADIR" -type f -size +0 | awk '{gsub(".*/","")}1' | awk '{gsub(".fasta","")}1' | sort)
## SGs in "$TADIR"
SG_PATTERNS=$(echo $PASSED_SG | awk '{gsub(" ","_|")}1' | awk '{gsub("super","")}1')
# Lines from "$OLD_SG" that correspond to SGs in $TADIR
SG_GROUPSTRINGS=$(grep -E "$SG_PATTERNS" <(cat "$GroupStringsPath_OLD" | sort) | awk '{gsub("_consensus","")}1' | sort -V)
# Save the 
echo "$SG_GROUPSTRINGS" > "$GroupStringsPath_NEW"
## remove $GroupStringsPath_OLD
rm "$GroupStringsPath_OLD"

