#!/bin/bash
#module load java
#module load python/3.7
# OldUsage: '/dipspades.sh' '/Path/To/Working/Directory/' '/Path/To/File_rename.csv' <nth individual to process>
# NewUsage: '/dipspades.sh' <nth individual to process> '/Path/To/Working/Directory/'
# Variable names for argument values
declare -i INDV=$1+1
WORKDIR=${2:-$PWD}

SAMPLEFILE=$WORKDIR'/File_rename.csv'
ASSEMBLYDIR=$WORKDIR'/Assembled_Contigs'

# Check that dipspades is the proper version, and if not exit with a message that version 12 must be used.
DIPSVER=$(dipspades.py --version)
DIPSEXP="SPAdes v3.12.0 [dipSPAdes mode]"
[ ! "$DIPSVER" = "$DIPSEXP" ] && echo "SPAdes v3.12.0 is required, older or newer versions will fail" && exit 1

# Path to dipspades.py
DIPSPADES='/panfs/pfs.local/home/j926w878/programs/SPAdes-3.12.0-Linux/bin/dipspades.py'
# DIPSPADES='dipspades.py'

# Getting the old and new patterns to use for filenames from a 2-column CSV file containing a header row
SAMPLEROW=$(head -n $INDV $SAMPLEFILE|tail -n 1)
OLDNAME=$(echo $SAMPLEROW | awk '{gsub(",.+", "")}1')
NEWNAME=$(echo $SAMPLEROW | awk '{gsub(".+,", "")}1')

# Where to process the samples
PROCDIR=$(echo $WORKDIR'/Processed_Samples')
SAMPLEDIR=$(echo $PROCDIR'/'$NEWNAME)

# If ASSEMBLYDIR does not exist, create it
[ ! -d "$ASSEMBLYDIR" ] && mkdir $ASSEMBLYDIR

# Move into ASSEMBLYDIR
cd $ASSEMBLYDIR

# Use dipspades to assemble preprocessed paired-end reads
$DIPSPADES --pe1-1 $SAMPLEDIR'/'$NEWNAME'_READ1_unmerged.fastq.gz' --pe1-2 $SAMPLEDIR'/'$NEWNAME'_READ2_unmerged.fastq.gz' --merged $SAMPLEDIR'/'$NEWNAME'_singletons.fastq.gz' -o $NEWNAME -k 21,33,55,77,99,127 --careful -t 10 -m 90 --hap-assembly --expect-gaps

echo $NEWNAME' dipspades Complete!'
