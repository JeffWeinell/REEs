#!/bin/bash
# usage: sbatch --nodes=1 --partition=sixhour --time=6:00:00 "/panfs/pfs.local/home/j926w878/work/conda/snakecap/test21_PhyloNet.sh" "/panfs/pfs.local/home/j926w878/programs/PhyloNet/InferNetwork_ML_bl_pl8_1_true.nex" "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/PhyloNet/test/InferNetwork_ML_bl_pl8_1_true_RUN2.out"
module load java
# module load anaconda
# conda activate py36

# Description: Runs ASTRAL III on IQTREE gene trees with low support nodes collapsed into polytomies
# 

# cp "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/gene-trees/supergroups_iqtree_OneSeqPerSample.trees" "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/supergroups_iqtree_OneSeqPerSample.trees"
# mkdir "/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/"

# Path to the PhyloNet jar file
PN="/panfs/pfs.local/home/j926w878/programs/PhyloNet/PhyloNet_3.8.2.jar"
# Path to the bin directory of Newick Utils
NU="/panfs/pfs.local/work/bi/bin/newick_utils/bin"
# Path to the input nexus file containing input data and parameters for the analysis
# INPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/ASTRAL/supergroups_iqtree_OneSeqPerSample.trees"
INPATH=$1
# Path where the output file should be saved
OUTPATH=$2
#OUTPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/PhyloNet/test.out"
MEM=${3:-"Xmx3000M"}

### Gene trees without treenames column
#TREES=$(awk '{print $2}' "$TREES_INPATH")
# Run ASTRAL on the uncollapsed gene trees, and allow java to use up to 3Gb of memory
# java -Xmx3000M -jar "$ASTRAL" -i <(echo "$TREES") -o "$TREE_OG_OUTPATH" 2>"$LOG_OG_OUTPATH"
# Graph result with Python as the root and save to 
# $NU/nw_reroot "$TREE_OG_OUTPATH" Python-molurus | $NU/nw_display -w 200 -S - > "$OUTDIR""/ASTRAL_1SeqPerSamplePerGene-BS0_graph.txt"

echo "Input: ""$INPATH"
echo "Output: ""$OUTPATH"

### Run PhyloNet
java -"$MEM" -jar "$PN" "$INPATH" > "$OUTPATH"
exit

#########################
#### PhyloNet Options ###
#########################
### Optional arguments are in square brackets
####
# InferNetwork_ML
####
# (geneTree1 [, geneTree2...]) The name of at least one gene tree is required, unless (all) is used meaning to use all gene trees.
# numReticulations The number of reticulations to include
# [-a taxa map] relates gene trees to species.
# [-bl] branch lengths should be considered
# [-b threshold]
# [-s startingNetwork]
# [-n numNetReturned]
# [-h {s1 [,s2...]}] 
# [-w (w1,w2,w3,w4)]
# [-x numRuns] 
# [-m maxNetExamined] 
# [-md moveDiameter] 
# [-rd reticulationDiameter]
# [-f maxFailure] 
# [-o] 
# [-po] 
# [-p (rel,abs)] 
# [-r maxRounds] 
# [-t maxTryPerBr] 
# [-i improveThreshold] 
# [-l maxBL] 
# [-pl numProcessors] 
# [-di] [result output file]
####
# MCMC_GT
####
# (geneTree1 [, geneTree2...]) The name of at least one gene tree is required, or (all) to use all gene trees.
# [-cl chainLength] 
# [-bl burnInLength] 
# [-sf sampleFrequency] 
# [-sd seed] 
# [-pp poissonParameter] 
# [-mr maximumReticulation] 
# [-pl parallelThreads] 
# [-tp temperatureList] 
# [-sn startingNetworkList] 
# [-tm taxonMap]
#

NEWICK_PATH=$1   # Character string with path to file containing newick trees (one per line).
PN_COMMANDS=$2 # Character string with PhyloNet commands.

NEWICK_PATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/PhyloNet/supergroups_iqtree_OneSeqPerSample_Resolved.trees"
PN_COMMANDS="InferNetwork_ML (all) 1 -bl -pl 8;"
NEXUS_PATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/PhyloNet/InferNetwork_ML_supergroups.nex"
OUTPATH="/panfs/pfs.local/home/j926w878/scratch/scratch_v4/SequenceCapture/PhyloNet/InferNetwork_ML_supergroups.out"

NEXUS=$(cat <(echo -e "#NEXUS\n\nBEGIN Trees;\n" ) <(awk '{print "Tree "$1"="$NF}' "$NEWICK_PATH") <(echo -e "\nEND;\n\nBEGIN PHYLONET; \n") <(echo -e "$PN_COMMANDS") <(echo -e "\nEND;"))
echo "$NEXUS" > "$NEXUS_PATH"

### Run PhyloNet
java -Xmx3000M -jar "$PN" "$NEXUS_PATH" > "$OUTPATH"

### Make a copy of the output network with inheritance probabilities removed, because Dendroscope and Icytree cannot read those.
# awk '{gsub("::[0-9.]+","")}1' "$OUTPATH" > "$OUTPATH""_NoInProbs"

# Constructing the NEXUS lines that PhyloNet reads.
#if [ $(awk '{if( NF > max ) max = NF} END {print max}' "$DATA_PATH") -eq 2 ]
#	then
#		NEXUS=$(cat <(echo -e "#NEXUS\n\nBEGIN Trees;\n" ) <(awk '{print "Tree "$1"="$NF}' "$DATA_PATH") <(echo -e "\nEND;\n\nBEGIN PHYLONET; \n") <(echo -e "$PN_COMMANDS") <(echo -e "\nEND;"))
#	#else 
#		# Just use the named newick trees...
#	#	exit 
#		# NEXUS=$(cat <(echo -e "#NEXUS\n\nBEGIN Trees;\n" ) <(awk '{print "Tree "$1"="$NF}' "$DATA_PATH") <(echo -e "\nEND;\n\nBEGIN PHYLONET; \n") <(echo -e "$PN_COMMANDS") <(echo -e "\nEND;"))
#fi

# echo 1:2:3:4:5 | awk -F: '{print $NF-1}'
# head -n 1
# echo -e "1:2:3:4:5\n1:2:3:4" | awk -F: '{print NF}' | sort -nu | tail -n 1
# echo -e "1:2:3:4:5\n1:2:3:4" | awk '{if( NF > max ) max = NF} END {print max}'









