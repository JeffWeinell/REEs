#!/bin/bash
module load 'blast+/2.9.0'
module load java
# Usage: '/targetMatching.sh' '/Path/To/Working/Directory/' '/Path/To/contigs.fa' '/Path/To/targetLoci.fa' '/Path/To/File_rename.csv' <nth individual to process>
# NewUsage: '/targetMatching.sh' '/Path/To/Working/Directory/' '/Path/To/targetLoci.fa' '/Path/To/File_rename.csv' <nth individual to process>
# NewUsage: '/targetMatching.sh' <nth individual to process> '/Path/To/Working/Directory/' '/Path/To/targetLoci.fa' '/Path/To/File_rename.csv'

# Variable names for argument values
declare -i INDV=$1+1
WORKDIR=${2:-'/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap2/'}
TARGETS=${3:-'/panfs/pfs.local/home/j926w878/work/SequenceCapture/Weinell_TargetLoci_Snakes_Final_18April2019.txt'}
SAMPLEFILE=${4:-'/panfs/pfs.local/home/j926w878/work/conda/snakecap/File_rename.csv'}

# WORKDIR='/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap2/'
# CONTIGS='/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap2/Processed_Samples/Psammodynastes-pulverulentus_KU328547/Psammodynastes-pulverulentus_KU328547_consensus-contigs-dipspades.fa'
# TARGETS='/panfs/pfs.local/home/j926w878/work/SequenceCapture/Weinell_TargetLoci_Snakes_Final_18April2019.txt'
# SAMPLEFILE='/panfs/pfs.local/home/j926w878/work/conda/snakecap/File_rename.csv'
# declare -i INDV=1+1

# Getting the sample name from the 2-column CSV file 'File_rename.csv'
SAMPLEROW=$(head -n $INDV $SAMPLEFILE|tail -n 1)
NEWNAME=$(echo $SAMPLEROW | awk '{gsub(".+,", "")}1')

# Where to process the samples
PROCDIR=$(echo $WORKDIR'/Processed_Samples')
SAMPLEDIR=$(echo $PROCDIR'/'$NEWNAME)

CONTIGS=$SAMPLEDIR'/'$NEWNAME'_consensus-contigs-dipspades.fa'

# Path to bbmap folder
BBMAP='/panfs/pfs.local/home/j926w878/programs/bbmap/'

# Make a blast database for the targets
## makeblastdb -in $TARGETS -parse_seqids -dbtype nucl -out $WORKDIR'/blast_db'

# DEDUPE the dipspades contigs...doesnt seem to have any affect except to make sequences interleaved.
cd $SAMPLEDIR
$BBMAP'dedupe.sh' in=$CONTIGS ordered=t overwrite=true out=$NEWNAME'_dd.fa' minidentity=97

blastn -task dc-megablast -db $WORKDIR'/blast_db' -query $SAMPLEDIR'/'$NEWNAME'_dd.fa' -out $SAMPLEDIR'/'$NEWNAME'_match.txt' -outfmt 6 -num_threads 10

echo $NEWNAME' target matching complete!'

