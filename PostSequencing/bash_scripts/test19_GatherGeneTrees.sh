#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test19_GatherGeneTrees.sh"
module load anaconda
conda activate py36

# Directory containing the treefiles, one per supergroup
TREESDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/gene-trees/SeqCap_TrimmedGenomes_OneSequencePerSample/treefiles"
# Path to the file 'allsamples_key.txt', which contains the sample name alias and original name in the first two columns, respectively.
SAMPLESKEYPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/allsamples_key.txt"
# Path where the output set of trees should be saved (two paths because two different naming schemes used for tips).
OUTPATH_TREES_1="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/gene-trees/supergroups_iqtree_OneSeqPerSample.trees"

## Path to the file 'allseqs_key.txt', which contains the sequence name alias and original name in the first two columns, respectively.
## SEQSKEYPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/allseqs_key.txt"
## OUTPATH_TREES_2="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/gene-trees/supergroups_genetrees_LongNames_iqtree.trees"

# Read all trees in $TREESDIR; for each tree, print the filepath and topology on separate lines, and then print a blank line; 
ALLTREES_0=$(find "$TREESDIR" -type f | xargs tail -n +1)
# Name to use for each tree, extracted from the the filepath
ALLTREES_NAMES=$(echo "$ALLTREES_0" | grep '.treefile' | awk '{print $2}' | awk '{gsub("^.*\\/","")}1' | awk '{gsub(".treefile","")}1' )
# Tree topologies (one per line)
ALLTREES_TOPOLOGY=$(echo "$ALLTREES_0" | awk NF | grep -v '.treefile')
# paste tree name and topology, respectively, separated a space
ALLTREES_1=$(paste -d " " <(echo "$ALLTREES_NAMES") <(echo "$ALLTREES_TOPOLOGY") | sort -V)

### Naming Scheme 1: Rename samples to their original names.
# remove sequence and group info 
ALLTREES_2=$(echo "$ALLTREES_1" | awk '{gsub("_R_","")}1' | awk '{gsub("_[0-9]*_group[0-9]*","_")}1' | awk '{gsub("_[0-9]*_supergroup[0-9]*","_")}1')
SAMPLESKEY=$(awk '{print $1"_", $2}' "$SAMPLESKEYPATH")
ALLTREES_3="$ALLTREES_2"
for X in $(seq 1 $(echo "$SAMPLESKEY" | wc -l)) ; do
	PATx=$(echo "$SAMPLESKEY" | awk -v x="$X" 'NR==x{print $1}')
	REPx=$(echo "$SAMPLESKEY" | awk -v x="$X" 'NR==x{print $2}')
	ALLTREES_3_temp=$(echo "$ALLTREES_3" | awk '{gsub("'$PATx'","'$REPx'")}1')
	ALLTREES_3="$ALLTREES_3_temp"
done

# Save the renamed set of trees with short names to $OUTPATH_TREES_1
echo "$ALLTREES_3" > "$OUTPATH_TREES_1"

########## Code below here needs to be updated before implementing this script on trees with multiple tips per sample

### Naming Scheme 2: Rename samples from the form "sample{N}_{SequenceID} to the form "Genus-species_VoucherID_ContigNumber_ContigWindow"
# remove sequence and group info 
## ALLTREES_2B=$(echo "$ALLTREES_1" | awk '{gsub("_R_","")}1' | awk '{gsub("_group[0-9]*","_")}1')
## ALLTREES_3B="$ALLTREES_2B"

## SEQSKEY=$(awk '{print $1"_", $2}' "$SEQSKEYPATH")

# Filter SEQSKEY to only include the names that exist in at least one of the gene trees...

## NAMES_0=$(echo "$SEQSKEY" | awk '{print $2}')
## USENAMES=$(echo "$SEQSKEY" | awk '{print $1}' | awk '{gsub("_contig_","\t")}1' | awk '{gsub("_length_","\t")}1' | awk '{gsub("_sliding:","\t")}1')
## DBNAMES2=$(cat <("$USENAMES"))
## echo "$NEWNAMES" | wc

## for X in $(seq 1 $(echo "$SEQSKEY" | wc -l)) ; do
## 	PATx=$(echo "$SEQSKEY" | awk -v x="$X" 'NR==x{print $1}')
## 	REPx=$(echo "$SEQSKEY" | awk -v x="$X" 'NR==x{print $2}')
## 	ALLTREES_3B_temp=$(echo "$ALLTREES_3B" | awk '{gsub("'$PATx'","'$REPx'")}1')
## 	ALLTREES_3B="$ALLTREES_3B_temp"
## done

## Save the renamed set of trees with long names (includes contig ID and window) to $OUTPATH_TREES_2
## echo "$ALLTREES_3B" > "$OUTPATH_TREES_2"

#TEST_1="$TEST"
#TEST_2B=$(echo "$TEST_1" | awk '{gsub("_R_","")}1' | awk '{gsub("_group[0-9]*","_")}1')



