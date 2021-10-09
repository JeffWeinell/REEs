#!/bin/bash
module load java
module load python/3.7
# Usage: '/preprocessingReads.sh' '/Path/To/Raw/Reads/Directory/' '/Path/To/Working/Directory/' '/Path/To/File_rename.csv' <nth individual to process>

# Variable names for argument values
RAWDIR=$1
WORKDIR=$2
SAMPLEFILE=$3
declare -i INDV=$4+1
#CONTAMINANTS=$5
#BBMAP=$6

# Define defaults for CONTAMINANTS and BBMAP. May implement this later.
# if [$ARG -eq '']; then echo 'success'
# fi

# Path to folder with contaminant genomes
CONTAMINANTS='/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap2/Contamination_Genomes'

# Path to bbmap folder
BBMAP='/panfs/pfs.local/home/j926w878/programs/bbmap/'

# RAWDIR='/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap2/raw_data'
# WORKDIR='/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap2/'
# SAMPLEFILE='/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap2/File_rename.csv'
# declare -i INDV=1+1

# Getting the old and new patterns to use for filenames from a 2-column CSV file containing a header row
SAMPLEROW=$(head -n $INDV $SAMPLEFILE|tail -n 1)
OLDNAME=$(echo $SAMPLEROW | awk '{gsub(",.+", "")}1')
NEWNAME=$(echo $SAMPLEROW | awk '{gsub(".+,", "")}1')

echo $NEWNAME' preprocessing initiated'

# Where to process the samples
PROCDIR=$(echo $WORKDIR'/Processed_Samples')
SAMPLEDIR=$(echo $PROCDIR'/'$NEWNAME)

# If WORKDIR does not exist, create it
[ ! -d "$WORKDIR" ] && mkdir $WORKDIR

# If PROCDIR does not exist, create it
[ ! -d "$PROCDIR" ] && mkdir $PROCDIR

# If SAMPLEDIR does not exist, create it
[ ! -d "$SAMPLEDIR" ] && mkdir $SAMPLEDIR

# If a copy of SAMPLEFILE does not already in WORKDIR, copy it there and give it the name 'File_rename.csv'
CPSAMPLEFILE=$WORKDIR'/File_rename.csv'
[ ! -f "$CPSAMPLEFILE" ] && cp $SAMPLEFILE $CPSAMPLEFILE

# Paths to original (not copied) raw reads
RAWREAD1=$RAWDIR'/'$(ls $RAWDIR | grep $OLDNAME.* | head -n 1)
RAWREAD2=$RAWDIR'/'$(ls $RAWDIR | grep $OLDNAME.* | tail -n 1)

# Paths where raw reads are to be copied
NEWPATH1=$SAMPLEDIR"/"$NEWNAME'_READ1.fastq.gz'
NEWPATH2=$SAMPLEDIR"/"$NEWNAME'_READ2.fastq.gz'

# Copying the raw reads to their new location
cp $RAWREAD1 $NEWPATH1
cp $RAWREAD2 $NEWPATH2

# Move into SAMPLEDIR
cd $SAMPLEDIR

# Use fastp to trim adaptors and filter low quality reads from copied raw reads. Currently fastp is in the system path when using snakecap conda environment
# A read is considered low quality if >40% of the bases have a phred quality score < Q15.
fastp --in1 $NEWNAME'_READ1.fastq.gz' --in2 $NEWNAME'_READ2.fastq.gz' --out1 $NEWNAME'_READ1_fastp.fastq.gz' --out2 $NEWNAME'_READ2_fastp.fastq.gz' --length_required 30 --low_complexity_filter --complexity_threshold 30 --html adapter_trimming_fastp.html --json adapter_trimming_fastp.json --report_title $NEWNAME --thread 10
echo $NEWNAME' trimming with fastp Complete!'

## give bbmap bash scripts executable privileges
chmod +x $BBMAP'/bbsplit.sh'
chmod +x $BBMAP'/dedupe.sh'
chmod +x $BBMAP'/reformat.sh'
chmod +x $BBMAP'/bbmerge-auto.sh'
chmod +x $BBMAP'/bbmerge.sh'

## Use bbsplit.sh from bbmap to remove contaminant reads
$BBMAP'/bbsplit.sh' -Xmx80g in1=$NEWNAME'_READ1_fastp.fastq.gz' in2=$NEWNAME'_READ2_fastp.fastq.gz' ref=$CONTAMINANTS outu1=$NEWNAME'_Cleaned_READ1.fastq.gz' outu2=$NEWNAME'_Cleaned_READ2.fastq.gz'

## dedupe duplicate removal. Only exact duplicates removed here.
$BBMAP'/dedupe.sh' -Xmx80g in1=$NEWNAME'_Cleaned_READ1.fastq.gz' in2=$NEWNAME'_Cleaned_READ2.fastq.gz' ordered=t overwrite=true out=$NEWNAME'_deduped_reads.fastq.gz' minidentity=100 ac=f

## reformat the deduped reads. READ1 and READ2 of each pair moved to separate files
$BBMAP'/reformat.sh' in=$NEWNAME'_deduped_reads.fastq.gz' out1=$NEWNAME'_deduped_READ1.fastq.gz' out2=$NEWNAME'_deduped_READ2.fastq.gz' int=t

####################
## paired-end read merging. Merges (mergeable) read pairs by overlap detection... Not sure why they would overlap...
# Identical to bbmerge.sh when -Xmx flag is used, but when the flag is not used bbmerge-auto.sh attempts to use all available memory, unlike bbmerge.sh
# Settings implemented #
# 'rem' Do not merge if the predicted insert size differs before and after extension.
# 'ecct': (='error-correct with Tadpole') algorith, which is more computational intensive but more accurate that the other algorithms.
# 'k=60': kmer overlap set to 60 (slower but more accurate)
# 'extend2=60': reads are extended by 60 only after a failed merge attempt
# No sequence trimming done at this step.
$BBMAP'/bbmerge-auto.sh' -Xmx80g in1=$NEWNAME'_deduped_READ1.fastq.gz' in2=$NEWNAME'_deduped_READ2.fastq.gz' verystrict=t rem k=60 extend2=60 ecct outu1=$NEWNAME'_READ1_unmerged.fastq.gz' outu2=$NEWNAME'_READ2_unmerged.fastq.gz' out=$NEWNAME'_merged_pe.fastq.gz'

## dedupe almost-exact duplicate removal of merged reads. Reads >=98% similar are removed.
# May want to change this, because some rare alleles might be lost? Need to look into this.
$BBMAP'/dedupe.sh' -Xmx80g in=$NEWNAME'_merged_pe.fastq.gz' ordered=t overwrite=true out=$NEWNAME'_singletons.fastq.gz' minidentity=98

### Write table with read removal stats. Only counts forward reads because the number of forward and reverse reads will be the same even though file sizes may be different.
# Initial number of forward and reverse reads. 
declare -i NREADSINITIAL=$(zcat < $NEWPATH1 | wc -l)/4
# Number of reads remaining after fastp adapter removal
declare -i NREADSFASTP=$(zcat < $SAMPLEDIR'/'$NEWNAME'_READ1_fastp.fastq.gz' | wc -l)/4
# Number of reads removed as a result of adapter trimming with fastp
declare -i NREMOVEDFASTP=$NREADSINITIAL-$NREADSFASTP
# Number of reads after removing contaminant reads
declare -i NREADSCLEAN=$(zcat < $SAMPLEDIR'/'$NEWNAME'_Cleaned_READ1.fastq.gz' | wc -l)/4
# Number of contaminant reads removed by bbsplit.sh
declare -i NREMOVEDCONTAMINANTS=$NREADSFASTP-$NREADSCLEAN
# Number of unique reads after exact duplicates removed
declare -i NREADSUNIQUE=$(zcat < $SAMPLEDIR'/'$NEWNAME'_deduped_READ1.fastq.gz' | wc -l)/4
# Number of exact duplicates removed by dedupe
declare -i NDUPSREMOVED=$NREADSCLEAN-$NREADSUNIQUE
# Number of read pairs merged with bbmerge
declare -i MERGED=$(zcat < $SAMPLEDIR'/'$NEWNAME'_merged_pe.fastq.gz' | wc -l)/4
# Number of unmerged read pairs after using bbmerge
declare -i UNMERGED=$(zcat < $SAMPLEDIR'/'$NEWNAME'_READ1_unmerged.fastq.gz' | wc -l)/4
# Number of 'singleton' merged reads after removing non-exact duplicates (>=98% identical)
declare -i NSINGLETONS=$(zcat < $SAMPLEDIR'/'$NEWNAME'_singletons.fastq.gz' | wc -l)/4
# Number of non-exact duplicates removed from the set of merged reads
declare -i NMERGEDREMOVED=$MERGED-$NSINGLETONS
# Number of reads remaining (unmerged pairs or merged singletons) 
declare -i NREADSFINAL=$UNMERGED+$NSINGLETONS
# Total reads removed during preprocessing
declare -i TOTALREMOVED=$NREMOVEDFASTP+$NREMOVEDCONTAMINANTS+$NDUPSREMOVED+$NMERGEDREMOVED

### write table with the read removal and read merging stats
echo "Source, nPairs, nPairsRemoved, nUnmerged, nMerged" > $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'
printf "initial raw reads, %s, 0, %s, 0\n" $NREADSINITIAL $NREADSINITIAL >> $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'
printf "fastp adapter-complexity filter, %s, %s, %s, 0\n" $NREADSFASTP $NREMOVEDFASTP $NREADSFASTP >> $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'
printf "bbsplit decontamination, %s, %s, %s, 0\n" $NREADSCLEAN $NREMOVEDCONTAMINANTS $NREADSCLEAN >> $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'
printf "dedupe exact duplicate removal, %s, %s, %s, 0\n" $NREADSUNIQUE $NDUPSREMOVED $NREADSUNIQUE >> $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'
printf "bbmerge read merging, %s, 0, %s, %s\n" $NREADSUNIQUE $UNMERGED $MERGED >> $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'
printf "dedupe non-exact duplicate removal, %s, %s, %s, %s\n" $NREADSFINAL $NMERGEDREMOVED $UNMERGED $NSINGLETONS >> $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'
printf "final counts and total removed, %s, %s, %s, %s\n" $NREADSFINAL $TOTALREMOVED $UNMERGED $NSINGLETONS >> $SAMPLEDIR'/'$NEWNAME'_preprocessing_stats.csv'

echo $NEWNAME' Complete!'

