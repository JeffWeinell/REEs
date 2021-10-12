#!/bin/bash
module load 'blast+/2.9.0'
# module load java
# Old Usage: '/targetMatching.sh' <nth individual to process> <'/Path/To/Working/Directory/'> <'/Path/To/DirectoryDB'> ['/Path/To/targetLoci.fa'] ['/Path/To/File_rename.csv']
# Current usage example: sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi  03.sh <'./blast_database3'> './ConsensusHaplocontigs_AllSamples_PlusTargets.fa' './ConsensusHaplocontigs_AllSamples_PlusTargets_matches_v2.txt' 
# Future Usage (not yet implemented): '/03.sh' <'/Path/To/BlastDB/Directory/'> <'/Path/To/seqFile'> ['options from blast3.sh']

BLASTDBDIR=$1         # Path to directory where the blast database is located. This script assumes that a directory only contains one database and that the database is named 'blast_db'
SEQFILE=$2
SEQNAMES=$(awk -F "\"* \"*" '{print $1}' $SEQFILE)
SEQPATHS=$(awk -F "\"* \"*" '{print $2}' $SEQFILE)
declare -i NUMSAMPLES=$(echo $SEQNAMES | wc -w)
# Create a directory to hold the blast directory if such a directory doesnt exist.
[ ! -d $BLASTDBDIR ] && mkdir $BLASTDBDIR
# IF a BLAST database does not already exist in BLASTDBDIR, then do the following:

if [ ! -f $BLASTDBDIR'/blast_db.nhr' ] || [ ! -f $BLASTDBDIR'/blast_db.nin' ] || [ ! -f $BLASTDBDIR'/blast_db.nog' ] || [ ! -f $BLASTDBDIR'/blast_db.nsd' ] || [ ! -f $BLASTDBDIR'/blast_db.nsi' ] || [ ! -f $BLASTDBDIR'/blast_db.nsq' ]; then
	echo "BLAST database does not exist. Creating one from fasta files listed in in seqFile."
	
	TEMPDIR=$BLASTDBDIR'/tempseqs'
	[ ! -d $TEMPDIR ] && mkdir $TEMPDIR
	### Check that all SEQPATHS exist and Generate an error if any doesnt exist
	# for i in $(seq 1 $NUMSAMPLES); do
	# 	SAMPLEi=$(echo $SEQNAMES | awk -v i="$i" '{print $i}')
	# 	INPATHi=$(echo $SEQPATHS | awk -v i="$i" '{print $i}')
	# 	# [ ! -f $INPATHi ] && echo $INPATHi' does not exist' && exit 1
	# done
	### Copy all sequence files listed in SEQFILE to a temporary directory within the directory where the blast database will be created
	for i in $(seq 1 $NUMSAMPLES); do
		SAMPLEi=$(echo $SEQNAMES | awk -v i="$i" '{print $i}')
		INPATHi=$(echo $SEQPATHS | awk -v i="$i" '{print $i}')
		TEMPPATHi=$TEMPDIR'/'$SAMPLEi'.fa'
		[ ! -f $TEMPPATHi ] && cp $INPATHi $TEMPPATHi
		### Line below is similar to copying $INPATHi to $TEMPPATHi, except only the first field (before whitespace) of sequence names is used.
		[ ! -f $TEMPPATHi ] && awk '{print $1}' $INPATHi > $TEMPPATHi
	done
	
	### Update the sequence names by appending the sample name to each sequence name of the ith sample
	for i in $(seq 36 $NUMSAMPLES); do
		# Text that will be appended to each sequence
		SAMPLEi=$(echo $SEQNAMES | awk -v i="$i" '{print $i}')
		# File that will be updated
		TEMPPATHi=$TEMPDIR'/'$SAMPLEi'.fa'
		# edits the name (NEED TO UPDATE THIS SO THAT NAMES ARE FIRST MODIFIED TO KEEP EVERYTHING BEFORE THE FIRST WHITESPACE)
		sed -i "/^>/s/$/_$SAMPLEi/" $TEMPPATHi
	done
	### concatenate all files listed in SEQFILE to 'allseqs.fa' in BLASTDBDIR
	# Paths to the copied and updated files
	TEMPPATHS=$(find $TEMPDIR -type f | sort)
	ALLSEQS=$BLASTDBDIR'/allseqs.fa'
	[ ! -f $ALLSEQS ] && cat $TEMPPATHS > $ALLSEQS
	# Make the blast database
	makeblastdb -in $ALLSEQS -parse_seqids -dbtype nucl -out $BLASTDBDIR'/blast_db'
else
	echo "Using existing database found in BLASTDBDIR"
fi


# Create a directory to hold a copy of the sequence files listed in the SEQFILE

# seqFile is a two column file with the names of samples in the first column and paths to fasta file in the second column. These sequences will be included in the blast database.
## seqFile="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/seqFile_SnakeSamplesAndGenomes_blast.txt"

## WORKDIR='/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples'
## PROCDIR=$WORKDIR'/Processed_Samples'
## TARGS=$WORKDIR'/Weinell_TargetLoci_Snakes_Final_18April2019.txt'
## # Where gathered sequences are/should be saved
## TEMPDIR=$WORKDIR'/tempseqs2'
## [ ! -d $TEMPDIR ] && mkdir $TEMPDIR
## # Where to save all seqs
## ALLSEQS=$WORKDIR'/ContigsMasked_AllSamples_PlusTargets.fa'

## cd $WORKDIR
## # Copy all '*_consensus-contigs-dipspades.fa' files into TEMPDIR
## SAMPLESEQS=$(find $PROCDIR -name "*_consensus-contigs-dipspades_masked-slow_min100bpValid.fa")
## # Number of samples
## declare -i NSAMP=$(echo $SAMPLESEQS | wc -w)
## for i in $(seq 1 $NSAMP); do 
##    SAMPLESEQSi=$(echo $SAMPLESEQS | awk -v i="$i" '{print $i}')
##    cp $SAMPLESEQSi $TEMPDIR
## done

## # append the sample name to each sequence name of the ith sample
## NEWSAMPLESEQS=$TEMPDIR'/'$(ls $TEMPDIR | sort)
## declare -i NSAMP=$(echo $NEWSAMPLESEQS | wc -w)
## for i in $(seq 1 $NSAMP); do 
##    NEWSAMPLESEQSi=$(echo $NEWSAMPLESEQS | awk -v i="$i" '{print $i}')
##    NAMEi=$(echo $(basename $(echo $NEWSAMPLESEQSi)) | awk -F'_' '{print $1 FS $2}')
##    # edits the name
##    sed -i "/^>/s/$/_$NAMEi/" $NEWSAMPLESEQSi
## done

## checking the first line to ensure that sequence names are correct
## for i in $(seq 1 $NSAMP); do 
##    NEWSAMPLESEQSi=$(echo $NEWSAMPLESEQS | awk -v i="$i" '{print $i}')
##    head -n 1 $NEWSAMPLESEQSi
## done
# head -n 1 "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database4/tempseqs/Hydrophis-curtus.fa"

### combine all of the sequences in TEMPDIR plus the targets into a single file
## cat $(find $TEMPDIR -type f) $TARGS > $ALLSEQS

### modifying sequence names further because the longest name must be 50 characters or less
## cat $ALLSEQS | awk '{gsub("_contig_.+_length_", "contig_")}1' > $ALLSEQS'_2'
## cat $ALLSEQS'_2' | awk '{gsub("-preocularis","")}1' > $ALLSEQS
## rm $ALLSEQS'_2'

### Create a blast directory in BLASTDBDIR to hold all sequences
# If BLASTDBDIR does not exist, create it and then make a blast database in the directory for the targets
## BLASTDBDIR="$1"             # Path to directory where the blast database is located. This script assumes that a directory only contains one database and that the database is named 'blast_db'
## [ ! -d "$BLASTDBDIR" ] && mkdir $BLASTDBDIR && makeblastdb -in $ALLSEQS -parse_seqids -dbtype nucl -out $BLASTDBDIR'/blast_db'

### Parameters passed as arguments. Some are optional and include defaults.
# ALLSEQS=$2            # Path the fasta file containing query sequences

###
MATCHES=$3            # Path where the output file (hit table in format 6) should be written.
EVAL="${4:-'10'}"     # Expected value (e-value). Default 10 (same as NCBI default). Decrease EVAL to return fewer but better hits, but setting this too low might cause some homologous regions to be missed.
MAXHSPs="${5:-'1'}"   # Maximum number of hits to return between a particular query-subject pair. NCBI default is 10, but here the default is 1.
WORDSIZE="${6:-'11'}" # Minimum number of contiguous, identically matching base pairs to return a hit. Contiguous bases can be separated by gaps when dc-megablast is used; word_size must be 11 or 12 for dc-megablast, ...
TASK="${7:-'dc-megablast'}" # The algorithm to use. Default is 'dc-megablast'. Other nucleotide options are 'blastn', 'megablast', 'tblastn', and 'tblastx'.
### Parameters not passed as arguments
GAPOPEN='5'         # Penalty multiplier for opening a gap. Increasing the value would likely return a set of hits with fewer gaps on avaerage, and possibly fewer total hits.
GAPEXTEND='2'       # Penalty multiplier for increasing the length of a gap. Increasing the value would likely return a set of hits with shorter gaps on avaerage, and possibly fewer total hits.
REWARD='2'          # reward value for each bp match
PENALTY='-3'        # penalty value for each bp mismatch
#DUST='20 64 1'      # Change value to '-dust' if you want to filter low-complexity query sequences. See publication here: https://pubmed.ncbi.nlm.nih.gov/16796549/
PERCIDENTITY='0'    # Percent identity cutoff. Minimum percent identity required between hits.
# WINDOWSIZE        # window_size; The minimum window size in which multiple hits can occur. NCBI default is 40, meaning that 39bp window cannot can contain two hits (between the same query and subject) but a 40bp window can.
blastn -task $TASK -db $BLASTDBDIR'/blast_db' -query $ALLSEQS -out $MATCHES -max_hsps $MAXHSPs -evalue $EVAL -outfmt 6 -num_threads 10 -word_size $WORDSIZE -gapopen $GAPOPEN -gapextend $GAPEXTEND -reward $REWARD -penalty $PENALTY -perc_identity $PERCIDENTITY

###### Modify this so that the consensus haplocontigs from each specimen or target sequences blast against the database containing the combined set of sample consensus contigs and targets.
###### Then, combine the output hit tables and create a network from the edge-matrix with quer name and subject name columns of the hit table. Align the contigs belonging to each distinct cluster of the network

# Variable names for argument values
# declare -i INDV=$1+1
# WORKDIR=${2:-'/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/'}
# TARGETS=${3:-$WORKDIR'/Weinell_TargetLoci_Snakes_Final_18April2019.txt'}
# SAMPLEFILE=${4:-$WORKDIR'/File_rename.csv'}
# 
# # Where to hold the BLAST database
# BLASTDBDIR=$WORKDIR'/blast_database2'
# 
# # Getting the sample name from the 2-column CSV file 'File_rename.csv'
# SAMPLEROW=$(head -n $INDV $SAMPLEFILE|tail -n 1)
# NEWNAME=$(echo $SAMPLEROW | awk '{gsub(".+,", "")}1')
# 
# # Where to process the samples
# PROCDIR=$(echo $WORKDIR'/BlastResults/')
# # SAMPLEDIR=$(echo $PROCDIR'/'$NEWNAME)
# SAMPLEDIR=$(echo $WORKDIR'/tempseqs/'$NEWNAME)
# 
# CONTIGS=$SAMPLEDIR'/'$NEWNAME'_consensus-contigs-dipspades.fa'
# 
# # Path to bbmap folder
# # BBMAP='/panfs/pfs.local/home/j926w878/programs/bbmap/'
# 
# # If BLASTDBDIR does not exist, create it and then make a blast database in the directory for the targets
# [ ! -d "$BLASTDBDIR" ] && mkdir $BLASTDBDIR && makeblastdb -in $TARGETS -parse_seqids -dbtype nucl -out $BLASTDBDIR'/blast_db'
# 
# # DEDUPE the dipspades contigs... this is just a hack to get sequences in an interleaved format that blast can use.
# # cd $SAMPLEDIR
# # $BBMAP'dedupe.sh' in=$CONTIGS ordered=t overwrite=true out=$NEWNAME'_dd.fa' minidentity=97
# 
# blastn -task dc-megablast -db $BLASTDBDIR'/blast_db' -query $SAMPLEDIR'/'$NEWNAME'_dd.fa' -out $SAMPLEDIR'/'$NEWNAME'_match.txt' -outfmt 6 -num_threads 10
# 
# echo $NEWNAME' target matching complete!'

# BLASTDBDIR="$1"           # Path to directory where the blast database is located. This script assumes that a directory only contains one database and that the database is named 'blast_db'
# QUERY="$2"                # Path to the fasta file containing query sequences
# MATCHES="$3"              # Path where the output file (hit table in format 6) should be written.
# EVAL="${4:-'10'}"         # Expected value (e-value). Default 10 (same as NCBI default). Decrease EVAL to return fewer but better hits, but setting this too low might cause some homologous regions to be missed.
# MAXHSPs="${5:-'1'}"       # Maximum number of hits to return between a particular query-subject pair. NCBI default is 10, but here the default is 1.
# WORDSIZE="${6:-'11'}"     # Minimum number of contiguous, identically matching base pairs to return a hit. Contiguous bases can be separated by gaps when dc-megablast is used; word_size must be 11 or 12 for dc-megablast, ...
# TASK="${7:-'dc-megablast'}" # The algorithm to use. Default is 'dc-megablast'. Other nucleotide options are 'blastn', 'megablast', 'tblastn', and 'tblastx'.
# ### Parameters not passed as arguments
# GAPOPEN='5'         # Penalty multiplier for opening a gap. Increasing the value would likely return a set of hits with fewer gaps on avaerage, and possibly fewer total hits.
# GAPEXTEND='2'       # Penalty multiplier for increasing the length of a gap. Increasing the value would likely return a set of hits with shorter gaps on avaerage, and possibly fewer total hits.
# REWARD='2'          # reward value for each bp match
# PENALTY='-3'        # penalty value for each bp mismatch
# #DUST='20 64 1'     # Change value to '-dust' if you want to filter low-complexity query sequences. See publication here: https://pubmed.ncbi.nlm.nih.gov/16796549/
# PERCIDENTITY='0'    # Percent identity cutoff. Minimum percent identity required between hits.
# # WINDOWSIZE        # window_size; The minimum window size in which multiple hits can occur. NCBI default is 40, meaning that 39bp window cannot can contain two hits (between the same query and subject) but a 40bp window can.
# blastn -task $TASK -db $BLASTDBDIR'/blast_db' -query $CONTIGS -out $MATCHES -max_hsps $MAXHSPs -evalue $EVAL -outfmt 6 -num_threads 10 -word_size $WORDSIZE -gapopen $GAPOPEN -gapextend $GAPEXTEND -reward $REWARD -penalty $PENALTY -perc_identity $PERCIDENTITY # -dust $DUST

###################### Not yet implemented. For now use the R script 'targetMatchingAssessment.R' after running this sh script on all samples.
# Sort contigs into 'unmatched', 'only weak matches', 'multiple moderate/strong matches', 'one strong match' categories
# subset of the table containing weak matches
# awk ' $4 < 40 || $11 >= 0.05 || $12 < 50 ' $SAMPLEDIR'/'$NEWNAME'_match.txt' > $SAMPLEDIR'/'$NEWNAME'_Hits-WeakMatches.tsv'
# code below here needs work.
# GOODMATCHES=$(awk ' $4 >= 40 && $11 < 0.05 && $12 >= 50 ' $WORKDIR'/dummy.tsv')
# awk $GOODMATCHES {$print $1} | uniq -c
# GOODMATCHES
# =$(awk ' $4 >= 40 && $11 < 0.05 && $12 >= 50 ' $WORKDIR'/dummy.tsv'
# TEST=$(cat $SAMPLEDIR'/'$NEWNAME'_match.txt')
# TEST2=$(head $SAMPLEDIR'/'$NEWNAME'_match.txt')
# head $(echo $TEST)
# awk -F"\t" '{print $1}' 
# TEST=$(awk -F"," '{print $1}' $SAMPLEFILE)
# awk ' $4 >= 40, $11 <=0.5 ' $WORKDIR'/dummy.tsv'

