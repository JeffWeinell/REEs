#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --mem=50Gb --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test16.sh"
# module load R
module load anaconda
conda activate py36

############
## Description: identifies which groups from /groups_fasta/ are not contained within supergroups that include //
## a genome-sampled sequence (i.e., groups not in "supergroupAlignments_SeqCapTargRefConsensus_GenomeSeqs_TrimmedToConsensus"), and copies those groups to //
## /supergroups_fasta/ and calls them "supergroups". The supergroups with genome-sampled sequences are also in the /supergroups_fasta/ directory. Therefore, //
## the supergroups in the directory supergroups_fasta correspond to homologs (some may be paralogs though?).
############

TADIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCapTargRefConsensus_GenomeSeqs_TrimmedToConsensus"
#STATSMATPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_ConsensusSCTR_vs_GenomesUntrimmed_STATS.txt"
SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test16.sh"
#STATSMATPATH_TA="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_ConsensusSCTR_vs_GenomesTrimmed_STATS.txt"
#STATSMAT_TA=$(seqkit stats $(find "$TADIR" -type f -name *fasta | sort -V ) -e -b )

# File with each line containing the names of groups in the supergroup; each name is separated by a the vertical bar symbol "|"
#GroupStringsPath="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/GroupsStrings_Consensus-vs-Genomes_ConsensusSeqs.txt"
GroupStringsPath="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupStrings_groups.txt"
# Directory containing sequences belonging to each group
GROUPSEQSDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupseqs/groups_fasta"
# Path to the file 'allsamples_key.txt', which contains the sample name alias and original name in the first two columns, respectively.
SAMPLESKEYPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/allsamples_key.txt"
# Output directory for sequences (SeqCap and Trimmed Genome seqs) of each supergroup
SG_SEQSDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/supergroups_fasta"

[ ! -d "$SG_SEQSDIR" ] && mkdir "$SG_SEQSDIR"

######################

STRINGS="$(cat "$GroupStringsPath")"
NUMSUPERGROUPS="$(echo "$STRINGS" | wc -l)"
#SAMPLESKEY=$(cat "$SAMPLESKEYPATH")

V1=$1
V2=$2
MODE=${3:-"A"}
#MODE=${3:-"B"}

if [ $MODE = "A" ] ; then
   # exit
   [ "$V2" -gt "$NUMSUPERGROUPS" ] && V2="$NUMSUPERGROUPS"
   for SGi in $(seq "$V1" "$V2") ; do
      # Groups in the ith supergroup
      GPSi="$(echo "$STRINGS" | awk -v i="$SGi" 'NR==i' | awk '{gsub("_consensus","")}1' |  awk '{gsub("\\|","\n")}1')"
      # Number of groups in the ith supergroup
      NUMGPSi=$(echo "$GPSi" | wc -l)
#      SAMPLENAME_TargRef=$(echo "$SAMPLESKEY" | grep 'targetRef' | awk '{print $1}')
      # Get all SeqCap sequences included in each group
#      SGi_SC=$(cat  $(echo "$GPSi" | awk -v dir="$GROUPSEQSDIR" '{print dir"/"$1".fasta"}') | seqkit grep -v -r -p "$SAMPLENAME_TargRef")
      SGi_SC=$(cat  $(echo "$GPSi" | awk -v dir="$GROUPSEQSDIR" '{print dir"/"$1".fasta"}'))
      # Name of the current supergroup
      SGi_NAME="super"$(echo "$GPSi" | sort -V | awk 'NR==1')
      # Get all trimmed genomic sequences if any exist
      if [ -f "$TADIR""/""$SGi_NAME"".fasta" ]
      then 
         SGi_GEN="$(cat "$TADIR""/""$SGi_NAME"".fasta")"
         SG_SEQS="$(cat <(echo "$SGi_SC") <(echo "$SGi_GEN"))"
      else 
         SG_SEQS="$SGi_SC"
      fi
      # remove gaps and write the dataset to $SG_SEQSDIR
      echo "$SG_SEQS" | seqkit seq --remove-gaps -w 0 > "$SG_SEQSDIR""/""$SGi_NAME"".fasta"
      # print stats
      seqkit stats "$SG_SEQSDIR""/""$SGi_NAME"".fasta" -b -e 
   done
fi

# Mode B copies each unassigned group fasta to $SG_SEQSDIR and calls it a supergroup
if [ $MODE = "B" ] ; then
    SG_NAMES="$(echo "$STRINGS" | awk '{gsub("_consensus","")}1' |  awk '{gsub("\\|","\n")}1')"
    # list of group fasta files not represented in supergroups
    UNASSIGNED="$(comm -13 <( echo "$SG_NAMES" | awk '{print $1".fasta"}' | sort ) <( ls "$GROUPSEQSDIR" | grep ^group* | sort ) | sort -V )"
    NUMUNASSIGNED=$(echo "$UNASSIGNED" | wc -l)
    [ "$V1" -gt "$NUMUNASSIGNED" ] && exit
    [ "$V2" -gt "$NUMUNASSIGNED" ] && V2="$NUMUNASSIGNED"
    for Gi in $(seq "$V1" "$V2"); do
      UNASSIGNED_Gi=$(echo "$UNASSIGNED" | awk -v i="$Gi" 'NR==i')
      #UNASSIGNED_Gi=$(echo $UNASSIGNED | awk -v i="$Gi" '{print $i}')
      cp $GROUPSEQSDIR"/"$UNASSIGNED_Gi $SG_SEQSDIR"/super"$UNASSIGNED_Gi
    done
fi


