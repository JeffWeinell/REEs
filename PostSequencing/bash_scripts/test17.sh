#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --mem=50Gb --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test17.sh"
module load R
module load anaconda
conda activate py36

####
## Runs MAFFT on each supergroup (SeqCap data + possibly trimmed genomes) containing at least 4 samples and at most one sequence per sample.
####

# Maximum number of sequences per sample in input dataset to generate alignment.
MAXSEQSPERSAMPLE=1
# Minimum number of samples per input dataset to generate an alignment.
MINSAMPLES=4

# Input directory with SeqCap and Trimmed Genome seqs of each supergroup
SG_SEQSDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/supergroups_fasta"
# Output directory for alignments of datasets from $SG_SEQSDIR
SG_ALNDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes"
# Path to table with name, number of seqs, and min, mean, max sequence length for each supergroup
STATSMATPATH_SG="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroups_unaligned_SeqCap_GenomesTrimmed_STATS.txt"
# Path to the file 'allsamples_key.txt', which contains the sample name alias and original name in the first two columns, respectively.
SAMPLESKEYPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/allsamples_key.txt"

[ ! -d "$SG_ALNDIR" ] && mkdir "$SG_ALNDIR"

# TADIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCapTargRefConsensus_GenomeSeqs_TrimmedToConsensus"
#STATSMATPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_ConsensusSCTR_vs_GenomesUntrimmed_STATS.txt"
SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test17.sh"
#STATSMATPATH_TA="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_ConsensusSCTR_vs_GenomesTrimmed_STATS.txt"
#STATSMAT_TA=$(seqkit stats $(find "$TADIR" -type f -name *fasta | sort -V ) -e -b )

# File with each line containing the names of groups in the supergroup; each name is separated by a the vertical bar symbol "|"
# GroupStringsPath="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/GroupsStrings_Consensus-vs-Genomes_ConsensusSeqs.txt"
# Directory containing sequences belonging to each group
# GROUPSEQSDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupseqs/groups_fasta"

######################

# STRINGS="$(cat "$GroupStringsPath")"
# NUMSUPERGROUPS="$(echo "$STRINGS" | wc -l)"

V1=$1
V2=$2
#MODE=${3:-"A"}
#MODE=${3:-"B"}

# Make a stats table for unaligned input datasets
if [ ! -f "$STATSMATPATH_SG" ]
then
   exit
   #STATSMAT_SG=$(seqkit stats $(find "$SG_SEQSDIR" -type f -name *fasta | sort -V ) -e -b )
   #echo "$STATSMAT_SG" > "$STATSMATPATH_SG"
fi

STATSMAT_SG=$(cat "$STATSMATPATH_SG" | awk 'NR>1')
NUMSUPERGROUPS=$(echo "$STATSMAT_SG" | wc -l)

[ "$V1" -gt "$NUMSUPERGROUPS" ] && exit
[ "$V2" -gt "$NUMSUPERGROUPS" ] && V2="$NUMSUPERGROUPS"

SAMPLESKEY=$(cat "$SAMPLESKEYPATH")
SAMPLENAME_TargRef=$(echo "$SAMPLESKEY" | grep 'targetRef' | awk '{print $1}')

for SGi in $(seq "$V1" "$V2") ; do
   # Groups in the ith supergroup
   SGi_NAME=$(echo "$STATSMAT_SG" | awk -v i="$SGi" 'NR==i' | awk '{print $1}')
   # Where to save output alignment
   OUTPATHi="$SG_ALNDIR""/""$SGi_NAME"
   # Skip if alignment already exists
   [ -f "$OUTPATHi" ] && continue
   # Path to ith unaligned supergroup
   SGi_PATH="$SG_SEQSDIR""/""$SGi_NAME"
   # read in the unaligned dataset, without the consensusSCTR or targetRef sequences
   SGi_SEQS=$(seqkit seq "$SGi_PATH" | seqkit grep -r -p sample | seqkit grep -v -r -p "$SAMPLENAME_TargRef")
   SGi_SEQNAMES=$(echo "$SGi_SEQS" | seqkit seq -n | awk '{gsub("^_R_","")}1')
   # Names of seqs
   # SGi_SEQNAMES=$(seqkit seq "$SGi_PATH" -n | awk '{gsub("^_R_","")}1')
   # Skip to next SGi if fewer sequences than $MINSAMPLES
   #NUMSEQSi=$(echo "$STATSMAT_SG" | awk -v i="$SGi" 'NR==i' | awk '{print $4}')
   NUMSEQSi=$(echo "$SGi_SEQNAMES" | wc -l)
   [ "$NUMSEQSi" -lt "$MINSAMPLES" ] && continue
  # Names of sample associated with each sequence
   SGi_SEQS_SAMPLENAMES=$(echo "$SGi_SEQNAMES" | awk '{gsub("_.+","")}1')
   # Number of samples
   NUMSAMPLESi=$(echo "$SGi_SEQS_SAMPLENAMES" | sort | uniq | wc -l )
   # Skip to next SGi if fewer samples than $MINSAMPLES
   [ "$NUMSAMPLESi" -lt "$MINSAMPLES" ] && continue
   # Number of sequences per sample
   NUMSEQSPERSAMPLE=$(Rscript -e "sample=commandArgs(trailingOnly=TRUE);print.data.frame(as.data.frame(table(sample)),row.names=F)" $SGi_SEQS_SAMPLENAMES)
   SEQSPERSAMPLE_MAX=$(echo "$NUMSEQSPERSAMPLE" | awk 'NR>1{print $2}' | sort | uniq | tail -n 1)
   # Skip to the next SGi if $SEQSPERSAMPLE_MAX > $MAXSEQSPERSAMPLE
   [ "$SEQSPERSAMPLE_MAX" -gt "$MAXSEQSPERSAMPLE" ] && continue
   # Print SGi and SGi_NAME to screen
   echo "SGi:""$SGi"" ""$SGi_NAME"
   # Run MAFFT and save alignment to $SG_ALNDIR
   #mafft --auto --adjustdirection --quiet <("$SGi_PATH") > "$OUTPATHi" &
   mafft --auto --adjustdirection --quiet <(echo "$SGi_SEQS") | seqkit seq -w 0 --upper-case  > "$OUTPATHi" &
done
wait

# FILESi=$(echo "$STATSMAT_SG" | awk -v x="$V1" -v y="$V2" 'NR>=x && NR<=y{print $1}')
# OUTPATHSi=$(echo "$FILESi" | awk -v dir="$SG_ALNDIR" '{print dir"/"$1}')
# generates an error if output file doesnt exist
# seqkit stats $OUTPATHSi -b -e




