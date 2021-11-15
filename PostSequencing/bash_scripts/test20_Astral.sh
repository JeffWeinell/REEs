#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test20_Astral.sh"
module load java
# module load anaconda
# conda activate py36

# Description: Runs ASTRAL III on IQTREE gene trees with low support nodes collapsed into polytomies
# 

# cp "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/gene-trees/supergroups_iqtree_OneSeqPerSample.trees" "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/supergroups_iqtree_OneSeqPerSample.trees"
# mkdir "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/"

# Path to the Astral III jar file
ASTRAL="/panfs/pfs.local/work/bi/bin/Astral_5.6.1/astral.5.6.1.jar"
# Path to the bin directory of Newick Utils
NU="/panfs/pfs.local/work/bi/bin/newick_utils/bin"
# Path to directory where output files should be saved
OUTDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/"
# Path to the input gene trees file
TREES_INPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/supergroups_iqtree_OneSeqPerSample.trees"
# Path where the collapsed gene trees should be written. Low support nodes from the input tree are collapsed into polytomies and written here.
TREES_BS10_PATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/supergroups_iqtree_OneSeqPerSample-BS10.trees"
# Path to the output tree for input trees uncollapsed
TREE_OG_OUTPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/ASTRAL_1SeqPerSamplePerGene-BS0.tree"
# Path to the output logfile for input trees uncollapsed
LOG_OG_OUTPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/ASTRAL_1SeqPerSamplePerGene-BS0.log"
# Path to the output tree for input trees with nodes collapsed
TREE_BS10_OUTPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/ASTRAL_1SeqPerSamplePerGene-BS10.tree"
# Path to the output logfile for input trees with nodes collapsed
LOG_BS10_OUTPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/ASTRAL_1SeqPerSamplePerGene-BS10.log"

### Gene trees without treenames column
TREES=$(awk '{print $2}' "$TREES_INPATH")
# Run ASTRAL on the uncollapsed gene trees, and allow java to use up to 3Gb of memory
# java -Xmx3000M -jar "$ASTRAL" -i <(echo "$TREES") -o "$TREE_OG_OUTPATH" 2>"$LOG_OG_OUTPATH"
# Graph result with Python as the root and save to 
# $NU/nw_reroot "$TREE_OG_OUTPATH" Python-molurus | $NU/nw_display -w 200 -S - > "$OUTDIR""/ASTRAL_1SeqPerSamplePerGene-BS0_graph.txt"

### Testing ASTRAL with only the first 5000 gene trees:
# TREES_5000=$(echo "$TREES" | awk 'NR<=5000{print}')
### Run ASTRAL on the first 5000 of the uncollapsed gene trees
# java -Xmx3000M -jar "$ASTRAL" -i <(echo "$TREES_5000") -o "$OUTDIR""/first_5000.tree" 2>"$OUTDIR""/first_5000.log"

### Gene trees with nodes having less than 10% UfBoot collapsed into polytomies.
$NU/nw_ed <(echo "$TREES") 'i & b<=10' o > $TREES_BS10_PATH
### Run ASTRAL on the collapsed gene trees
java -Xmx3000M -jar "$ASTRAL" -i "$TREES_BS10_PATH" -o "$TREE_BS10_OUTPATH" 2>"$LOG_BS10_OUTPATH"
# Graph result with Python as the root and save to 
$NU/nw_reroot "$TREE_BS10_OUTPATH" Python-molurus | $NU/nw_display -w 200 -S - > "$OUTDIR""/ASTRAL_1SeqPerSamplePerGene-BS10_graph.txt"

## ### Score the output tree by comparing it to the set of input gene trees
## # For the fully resolved gene trees
## java -jar "$ASTRAL" -q "$TREE_OG_OUTPATH" -i <(echo "$TREES") -o "$OUTDIR""/OG_simulated_scored.tre" 2> "$OUTDIR""/OG_simulated_scored.log"
## # For the collapsed gene trees
## java -jar "$ASTRAL" -q "$TREE_BS10_OUTPATH" -i "$TREES_BS10_PATH" -o "$OUTDIR""/BS10_simulated_scored.tre" 2> "$OUTDIR""/BS10_simulated_scored.log"

######### Print the tree, rooting with Python
### $NU/nw_reroot "$OUTDIR""/first_500.tree" Python-molurus | $NU/nw_display -S -










