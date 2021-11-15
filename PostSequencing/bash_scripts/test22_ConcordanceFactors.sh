#!/bin/bash
# sbatch --nodes=1 --partition=sixhour --ntasks-per-node=8 --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test22_ConcordanceFactors.sh"
# module load java
# module load perl
module load anaconda
conda activate py36

# Description: Runs iqtree to calculate concordance factors given a set of gene trees and a reference tree.
# For now I am using the ASTRAL tree as the reference tree.
# Used concat.sh to concatenate gene alignments.

## REFTREE_PATH=$1
## GENETREES_PATH=$2 # Each line must start with a '('
## ALIGNMENT_PATH=$3
## OUTDIR=$4

REFTREE_PATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/ASTRAL_1SeqPerSamplePerGene-UFBoot0.tree"
GENETREES_PATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/supergroups_iqtree_OneSeqPerSample_Resolved.trees"
OUTDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/"
ALIGNMENT_PATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes_OGNames_concatenated.fasta"

# work around because my gene trees file has tree names in the first column
GENETREES=$(awk '{print $2}' "$GENETREES_PATH")

#iqtree -t "$REFTREE_PATH" --gcf "$GENETREES_PATH" -s "$ALIGNMENT_PATH" --scf 100 --prefix "$OUTDIR""/concord"
iqtree -t "$REFTREE_PATH" --gcf <(echo "$GENETREES") -s "$ALIGNMENT_PATH" --scf 100 --prefix "$OUTDIR""/concord"

# iqtree -t concat.treefile --gcf loci.treefile -s alignment.nex --scf 100 --prefix concord

##### Switching back to the old names
# INDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes"
# SAMPLESKEYPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/SnakeCap_AllSamples/blast_database10/allsamples_key.txt"
# OUTDIR1="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes_SampleName_ContigName"
# OUTDIR2="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes_SampleName"
# OUTDIR3="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCap_TrimmedGenomes_OGNames"
# 
# ALIGNMENT_PATHS=$(find "$INDIR" -type f -name "*.fasta")
# NUMPATHS=$(echo "$ALIGNMENT_PATHS" | wc -l)
# 
# OUTPATHS1=$(echo "$ALIGNMENT_PATHS" | awk '{gsub(".*/","")}1' | awk -v dir="$OUTDIR1" '{print dir"/"$1}')
# OUTPATHS2=$(echo "$ALIGNMENT_PATHS" | awk '{gsub(".*/","")}1' | awk -v dir="$OUTDIR2" '{print dir"/"$1}')
# OUTPATHS3=$(echo "$ALIGNMENT_PATHS" | awk '{gsub(".*/","")}1' | awk -v dir="$OUTDIR3" '{print dir"/"$1}')
# 
# ALIAS=$(awk 'NF=2' "$SAMPLESKEYPATH" | awk '{gsub(" ","\t")}1')
# 
# [ ! -d "$OUTDIR1" ] && mkdir "$OUTDIR1"
# [ ! -d "$OUTDIR2" ] && mkdir "$OUTDIR2"
# [ ! -d "$OUTDIR3" ] && mkdir "$OUTDIR3"
# 
# for X in $(seq 1 "$NUMPATHS") ; do 
#   INPATHi=$(echo "$ALIGNMENT_PATHS" | awk -v x="$X" 'NR==x')
#   OUTPATH_1i=$(echo "$OUTPATHS1" | awk -v x="$X" 'NR==x')
#   cat "$INPATHi" | awk '{gsub("_R_","")}1' > "$OUTPATH_1i" &
# done
# wait
# 
# for X in $(seq 1 "$NUMPATHS") ; do 
#   OUTPATH_1i=$(echo "$OUTPATHS1" | awk -v x="$X" 'NR==x')
#   OUTPATH_2i=$(echo "$OUTPATHS2" | awk -v x="$X" 'NR==x')
#   cat "$OUTPATH_1i" | awk '{gsub("_.*","")}1' > "$OUTPATH_2i" &
# done
# wait 
# 
# for X in $(seq 1 "$NUMPATHS") ; do 
#   OUTPATH_2i=$(echo "$OUTPATHS2" | awk -v x="$X" 'NR==x')
#   OUTPATH_3i=$(echo "$OUTPATHS3" | awk -v x="$X" 'NR==x')
#   seqkit replace -p '(.+)$' -r '{kv}' -k <(echo "$ALIAS") "$OUTPATH_2i" > "$OUTPATH_3i" &
# done
# wait
# 

