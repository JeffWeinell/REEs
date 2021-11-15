#!/bin/bash
############################################################
# SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test14.sh"
# INDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/supergroups_fasta_SeqCapTargRefConsensus_GenomesUntrimmed"
# OUTDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupAlignments_SeqCapTargRefConsensus_GenomesUntrimmed_max100seqs/"
# GROUPSFILE="NA"
#
# usage: sbatch --nodes=1 --partition=sixhour --time=6:00:00 "$SHij" "$INDIR" "$OUTDIR" "NA" "A" "NA" "NA"
#
### Description:
# Runs MAFFT on each FASTA file in $INDIR and write 
# Each input FASTA corresponds to a supergroup, and includes ConsensusSeqCapTargRef sequences and their untrimmed (20Kb windowed) contig matches of genomic samples.
# The alignments generated from this step will be trimmed to the span of the ConsensusSeqCapTargRef seqs; the trimmed sequences of genome samples will then be aligned with the SeqCap and TargRef sequences.
### NOTE: run this if you cause an infinite loop:  squeue -u $USER -h | awk '{print $1}' | xargs scancel
# module load R
module load anaconda
conda activate py36
INDIR=$1               # directory containing fasta files with unaligned sequences 
OUTDIR=$2              # directory where aligned sequences should be saved
# GROUPSFILE=${3:-"NA"}  # Optional file indicating which groups in $INDIR to evaluate
# MODE=${4:-"A"}
MODE=${3:-"A"}
SHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/test14.sh"
PAR=${5:-5}
MAXSEQS=100  # Alignment is skipped for input fasta files containing more sequences than MAXSEQS
# STARTj=${5:-"1"}

#[ "$GROUPSFILE" = "NA" ] && INPATHS=$(find $INDIR -type f -name "supergroup*")
INPATHS=$(find $INDIR -type f -name "supergroup*" | sort -V)
#[ ! "$GROUPSFILE" = "NA" ] && INPATHS=$(awk '{ print "'$INDIR'"$1".fasta"}' $GROUPSFILE)

NUMPATHS=$(echo $INPATHS | wc -w)

[ ! -d "$OUTDIR" ] && mkdir -p "$OUTDIR"

if [ $MODE = "A" ] ; then
  for VARj in $(seq 1 $PAR $NUMPATHS) ; do
  #for VARj in $(seq 1 $PAR 4) ; do 
    #sbatch --nodes=1 --ntasks-per-node=10 --mem=50Gb --time=6:00:00 --partition=sixhour $SHij $INDIR $OUTDIR $GROUPSFILE "B" $VARj
    sbatch --nodes=1 --ntasks-per-node=10 --mem=50Gb --time=6:00:00 --partition=sixhour $SHij $INDIR $OUTDIR "C" $VARj
    # sbatch --nodes=1 --partition=sixhour $SHij $INDIR $OUTDIR $GROUPSFILE "B" $VARj
  done
fi


if [ $MODE = "B" ] ; then
  VARj=$4
  VARj2=$(($VARj+$PAR-1))
  [ $VARj2 -ge $NUMPATHS ] && VARj2=$NUMPATHS
  #INPATHSi=$(echo "$INPATHS" | cut $"-f""$VARj""-""$VARj2" -d" ")
  INPATHSi=$(echo "$INPATHS" | awk -v varj="$VARj" -v varj2="$VARj2" 'NR>=varj && NR<=varj2 {print}')
  NUMPATHSi=$(echo "$INPATHSi" | wc -l)
  INFILESi=$(echo "$INPATHSi" | awk '{gsub(".+/","")}1')
  OUTPATHSi=$(echo "$INFILESi" | awk -v outdir="$OUTDIR" '{print outdir"/"$1}')
  INFOi=$(paste <(seq $VARj $VARj2) <(echo "$INFILESi"))
  # Print some info to screen.
  echo -e "VARi\tDataset\n""$INFOi""\n\nQuietly running MAFFT in parallel...\n\n"
  for VARi in $(seq 1 $NUMPATHSi) ; do
    INPATHi=$(echo "$INPATHSi" | awk -v i=$VARi 'NR==i{print}')
    NUMSEQSi=$(seqkit seq --name "$INPATHi" | wc -l)
    [ $NUMSEQSi -gt $MAXSEQS ] && continue
    OUTPATHi=$OUTDIR"/"$(echo $INPATHi | awk '{gsub(".+/","")}1')
    mafft --localpair --adjustdirection --quiet --maxiterate 1000 "$INPATHi" > "$OUTPATHi""_temp" &
  done
  wait
  
  for VARi in $(seq 1 $NUMPATHSi) ; do
    INPATHi=$(echo $INPATHSi | awk -v i=$VARi '{print $i}')
    OUTPATHi=$OUTDIR"/"$(echo $INPATHi | awk '{gsub(".+/","")}1')
    #OUTPATHi=$OUTDIR"/"$(basename $INPATHi)
    NUMSEQSi=$(seqkit seq --name "$INPATHi" | wc -l)
    [ $NUMSEQSi -gt $MAXSEQS ] && continue
    seqkit seq --upper-case --line-width 0 "$OUTPATHi""_temp" > "$OUTPATHi"
    rm "$OUTPATHi""_temp" &
  done
  wait
  #echo "Success."
fi


if [ $MODE = "C" ] ; then
  VARj=$4
  PAR=$5
  VARj2=$(($VARj+$PAR-1))
  [ $VARj2 -ge $NUMPATHS ] && VARj2=$NUMPATHS
  #INPATHSi=$(echo "$INPATHS" | cut $"-f""$VARj""-""$VARj2" -d" ")
  INPATHSi=$(echo "$INPATHS" | awk -v varj="$VARj" -v varj2="$VARj2" 'NR>=varj && NR<=varj2 {print}')
  NUMPATHSi=$(echo "$INPATHSi" | wc -l)
  INFILESi=$(echo "$INPATHSi" | awk '{gsub(".+/","")}1')
  OUTPATHSi=$(echo "$INFILESi" | awk -v outdir="$OUTDIR" '{print outdir"/"$1}')
  INFOi=$(paste <(seq $VARj $VARj2) <(echo "$INFILESi"))
  # Print some info to screen.
  echo -e "\nVARi\tDataset\n""$INFOi"
  for VARi in $(seq 1 $NUMPATHSi) ; do
    INPATHi=$(echo "$INPATHSi" | awk -v i=$VARi 'NR==i{print}')
    NUMSEQSi=$(seqkit seq --name "$INPATHi" | wc -l)
   # [ $NUMSEQSi -gt $MAXSEQS ] && continue
    OUTPATHi=$OUTDIR"/"$(echo $INPATHi | awk '{gsub(".+/","")}1')
    #SGi=$(basename "$INPATHi")
    #CONGENi_INPATH="$INDIR""/""$SGi"
#    SEQS_GEN_U=$(cat "$INPATHi" | seqkit grep --use-regexp --pattern 'sample')
#    SEQS_CON_U=$(cat "$INPATHi" | seqkit grep --use-regexp --pattern 'consensus')
#    ### Chop up long genome sequences into more manageable sizes and use a sliding window to be sure than fragments can find each other
#    SEQS_GEN_WIND_U=$(echo "$SEQS_GEN_U" | seqkit sliding --step 1800 --window 2000 --greedy | seqkit seq --min-len 100 --remove-gaps)
#    SEQS_CONU_GENWU=$(cat <(echo "$SEQS_CON_U") <(echo "$SEQS_GEN_WIND_U"))
#    #mafft --localpair --adjustdirection --quiet --maxiterate 1000 <(echo "$SEQS_CONU_GENWU") > "$OUTPATHi""_temp"  &
#    echo "$SEQS_CONU_GENWU" > "$OUTPATHi""_windows_temp"
#    #echo "$SEQS_CONU_GENWU" | seqkit stats
    echo -e "\nQuietly running MAFFT on dataset..."
#    seqkit stats $OUTPATHi""_windows_temp -b -e
    seqkit stats "$INPATHi" -b -e
    #mafft --localpair --adjustdirection --quiet --maxiterate 1000 "$OUTPATHi""_windows_temp" > "$OUTPATHi""_temp"  &
    #mafft --auto --adjustdirection --quiet "$OUTPATHi""_windows_temp" > "$OUTPATHi""_temp"  &
    mafft --auto --adjustdirection --quiet "$INPATHi" > "$OUTPATHi""_temp"  &
  done
  wait
  
  echo "alignment stats:"
  #OUTPATHSi_STRING=$(echo -n $(echo "$OUTPATHSi" | awk '{print $1"_temp"}') | awk '{gsub(" ",", ")}1')
  # OUTPATHSi_TEMP=$(echo "$OUTPATHSi" | awk '{print $1"_temp"}' | awk '{gsub("$"," ")}1')
  # OUTPATHSi_TEMP="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/test/supergroup1.fasta_temp /panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/test/supergroup3.fasta_temp"
  # Rscript -e "args<-commandArgs(trailingOnly=T); var=unlist(strsplit(args[1],split=\" \")); if any(file.exists(var)) { var[file.exists(var)]} " "$OUTPATHSi_TEMP"
  #OUTPATHSi_TEMP=$(echo "$OUTPATHSi" | awk '{print $1"_temp"}')
  STATPATHS1="$(echo "$OUTPATHSi" | awk '{print $1"_temp"}')"
  #seqkit stats "$OUTPATHSi"* -e -b
  seqkit stats $STATPATHS1 -e -b
# '{printf "%s+",$0} END {print ""}'
  # [ ! "$CHECK2" -lt 2 ] && continue
  for VARi in $(seq 1 $NUMPATHSi) ; do
    INPATHi=$(echo $INPATHSi | awk -v i=$VARi '{print $i}')
    OUTPATHi=$OUTDIR"/"$(echo $INPATHi | awk '{gsub(".+/","")}1')
#    NUMSEQSi=$(seqkit seq -n "$INPATHi" | wc -l)
#    [ $NUMSEQSi -gt $MAXSEQS ] && echo "More than ""$MAXSEQS" continue
    # Check that the output alignment exists and if not echo info and continue to next file
    [ ! -f "$OUTPATHi""_temp" ] && echo $(basename "$OUTPATHi"" failed") && continue
    # DONT USE REMOVE GAPS FLAG AFTER THE ALIGNMENT HAS BEEN GENERATED
    cat "$OUTPATHi""_temp" | seqkit seq --upper-case --line-width 0  > "$OUTPATHi" &
#    rm "$OUTPATHi""_temp" # &
  done
  wait
  # echo "Success."
  echo -e "\nfinal alignment stats:"
  STATPATHS2="$(cat <(echo "$OUTPATHSi") <(echo "$OUTPATHSi" | awk '{print $1"_temp"}'))"
  # seqkit stats "$OUTPATHSi"* -b -e
  # Dont put quotes around STATPATHS2
  seqkit stats $STATPATHS2 -b -e
  echo -e "\n"
fi

exit 

### Code below should be used in test16.sh to merge windows belonging to the same sample sequence
# .libPaths("/panfs/pfs.local/home/j926w878/work/R-packages")
# library(Biostrings)
# library(gtools)
# 
# args <- commandArgs(trailingOnly=TRUE)
# indir    <- args[1]
# outdir   <- args[2]
# varj     <- as.numeric(args[3])
# varj2    <- as.numeric(args[4])
# 
# print(data.frame(indir=indir,outdir=outdir,varj=varj,varj2=varj2))
# 
# infiles  <- gtools::mixedsort(list.files(indir,full.names=T))[varj:varj2]
# outfiles <- file.path(outdir,basename(infiles))
# 
# print("A")
# conNames      <- paste0(tools::file_path_sans_ext(basename(infiles)),"_consensus")
# print("B")
# aln.list      <- lapply(infile, Biostrings::readDNAMultipleAlignment)
# 
# ### Get set of sliding windows for a particular sample sequence
# print("C")
# ### Replace gaps with Ns
# aln2.list     <- lapply(aln.list,function(x){Biostrings::DNAMultipleAlignment(gsub("-","N",x))})
# print("D")
# conStr.list   <- lapply(aln2.list,function(x){Biostrings::consensusString(x, ambiguityMap=IUPAC_CODE_MAP)})
# print("E")
# conDNA.list   <- lapply(conStr.list,function(x){Biostrings::DNAStringSet(x)})
# print("F")
# conDNA        <- do.call(c,conDNA.list)
# print("G")
# names(conDNA) <- conNames
# print("F")
# out           <- lapply(1:length(outfiles),function(x){Biostrings::writeXStringSet(conDNA[x],outfiles[x])})
# 
# # Merge sliding window sequences derived from the same sample and contig. This is done by gettin the consensus.
# Biostrings::readDNAMultipleAlignment
# 
#INPATHi="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/supergroupseqs/supergroups_fasta_SeqCapTargRefConsensus_GenomesUntrimmed/supergroup5823.fasta"
#SEQS_CON_U=$( cat "$INPATHi" | seqkit grep -r -p 'consensus')
## Number of consensus seqs in the input dataset
#NSEQS_CON_U=$(echo $SEQS_CON_U | seqkit stats | awk 'NR==2{print $4}')
## If $NSEQS_CON_U = 1, then align the other sequences to $CON_U with the --addlong flag
#mafft --localpair --adjustdirection --quiet --maxiterate 1000 <(echo $SEQS_CON_U) > "$OUTPATHi""_temp1"
## If $NSEQS_CON_U > 1, then align the the $CON_U sequences, and then align the other sequences to the $ALN_CON
#mafft --localpair --adjustdirection --quiet --maxiterate 1000 <(echo $SEQS_CON_U) > "$OUTPATHi""_temp1"
#ALN_CON=$(cat "$OUTPATHi""_temp1")

# MEMJOBS=$(sacct -u $USER | awk 'NR>2{if(($1>=27832993) && ($2=="test14.sh") && ($6=="OUT_OF_ME+") ){print $1}}')


