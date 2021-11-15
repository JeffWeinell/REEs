#!/bin/bash
module load 'blast+/2.9.0'
module load anaconda
conda activate py36

# Defaults shown in curly brackets
# usage: sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi  "/panfs/pfs.local/home/j926w878/work/conda/snakecap/03-temp.sh" </path/to/blastDB/directory> </path/to/input/seqFile> </path/output/hitTable> [eval {10}] [max hits per query-subject pair {1}] [wordsize {11}] [task {dc-megablast}]
# sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi "/panfs/pfs.local/home/j926w878/work/conda/snakecap/03-temp.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database7" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/seqFile_test_blast.txt"
# sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi "/panfs/pfs.local/home/j926w878/work/conda/snakecap/03-temp.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database8" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/seqFile_SnakeSamplesAndGenomes_blast.txt"
# sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi "/panfs/pfs.local/home/j926w878/work/conda/snakecap/03-temp.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database9" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/seqFile_SnakeSamplesAndGenomes_blast.txt"
# sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi "/panfs/pfs.local/home/j926w878/work/conda/snakecap/03-temp.sh" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/seqFile_SnakeSamplesAndGenomes_blast.txt"

BLASTDBDIR=$1 # Path to directory where the blast database is located. This script assumes that a directory only contains one database and that the database is named 'blast_db'
SEQFILE=$2

# BLASTDBDIR="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10"
# SEQFILE="/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/seqFile_SnakeSamplesAndGenomes_blast.txt"

SEQNAMES=$(awk -F "\"* \"*" '{print $1}' $SEQFILE)
SEQPATHS=$(awk -F "\"* \"*" '{print $2}' $SEQFILE)
SEQTYPES=$(awk -F "\"* \"*" '{print $3}' $SEQFILE) # one of "SeqCap", "genome", or "targetRef", which refer to the source of the sequences.

declare -i NUMSAMPLES=$(echo $SEQNAMES | wc -w)
# Create a directory to hold the blast directory if such a directory doesnt exist.
[ ! -d $BLASTDBDIR ] && mkdir $BLASTDBDIR

## ##  # In the case that a blast database does not already exist, should a separate blast database be created for each sample and each sample blasted against every database except its own? Currently this must be false. When false, a single database is constructed (unless it already exists) from all sequences of all samples, and then all sequences of all samples are blasted against this database; BLAST may split up the database for some parallelization, but the better parallelization might be accomplished by setting PAIRWISE to true once this feature becomes available.
## ##  # PAIRWISE=false
## ##  # PAIRWISE=true

TEMPDIR=$BLASTDBDIR'/tempseqs'
[ ! -d $TEMPDIR ] && mkdir $TEMPDIR
# A directory hold intermediate copies of sequences.
TEMPDIR0=$TEMPDIR'0'
[ ! -d $TEMPDIR0 ] && mkdir $TEMPDIR0

### Copy all sequence files listed in SEQFILE to a temporary directory in the directory where the blast database will be created
# declare -i NUMSAMPLES=3
for i in $(seq 1 $NUMSAMPLES); do
# for i in $(seq 1 3); do
  SAMPLEi=$(echo $SEQNAMES | awk -v i=$i '{print $i}')
  ## echo $SAMPLEi
  INPATHi=$(echo $SEQPATHS | awk -v i=$i '{print $i}')
  TEMPPATHi=$TEMPDIR'/'$SAMPLEi'.fa'
  TEMPPATHi_0=$TEMPDIR0'/'$SAMPLEi'.fa'
  # Make a copy of the input sequences; while copying, break up any sequences >20Kb into windows.
  cat $INPATHi | seqkit sliding -W 20000 -s 20000 -g > $TEMPPATHi_0
  # rename sequences and do not wrap sequences.
  # cat $TEMPPATHi_0 | seqkit replace -p .+ -r 'sample'$i'_'{nr} --nr-width 10 > $TEMPPATHi
  cat $TEMPPATHi_0 | seqkit replace -p .+ -r 'sample'$i'_'{nr} --nr-width 10 -w 0 > $TEMPPATHi
done

### Modify the sequence names to have the form Sample{N}_{NR}, where {N} is the Nth sample and {NR} is the jth sequence of the Nth sample. Create a Key file to translate between oldnames and newnames.
KEY=$BLASTDBDIR'/allseqs_key.txt'
OLDNAMES=$BLASTDBDIR'/oldnames.txt'
NEWNAMES=$BLASTDBDIR'/newnames.txt'
SAMPLENAMES=$BLASTDBDIR'/samplenames.txt'
[ -f $OLDNAMES ] && rm $OLDNAMES
[ -f $NEWNAMES ] && rm $NEWNAMES
[ -f $SAMPLENAMES ] && rm $SAMPLENAMES
[ -f $KEY ] && rm $KEY
touch $OLDNAMES $NEWNAMES $SAMPLENAMES $KEY

#SAMPLE1=$(echo $SEQNAMES | awk -v i=1 '{print $i}')
#INPATH1=$(echo $SEQPATHS | awk -v i=1 '{print $i}')
#TEMPPATH1=$TEMPDIR'/'$SAMPLE1'.fa'
#TEMPPATH1_0=$TEMPDIR0'/'$SAMPLE1'.fa'
#echo $INPATH1
#echo $TEMPPATH1_0
#grep "^>" $TEMPPATH1_0 | awk '{gsub("^>","")}1' > $OLDNAMES

for i in $(seq 1 $NUMSAMPLES); do
  SAMPLEi=$(echo $SEQNAMES | awk -v i=$i '{print $i}')
  TEMPPATHi=$TEMPDIR'/'$SAMPLEi'.fa'
  TEMPPATHi_0=$TEMPDIR0'/'$SAMPLEi'.fa'
  # append the old sequence names to OLDNAMES
  #[ $i = 1 ] && grep "^>" $TEMPPATHi_0 | awk '{gsub("^>","")}1' > $OLDNAMES
  grep "^>" $TEMPPATHi_0 | awk '{gsub("^>","")}1' >> $OLDNAMES
  # Number of sequences of ith sample
  declare -i NUMSEQSi=$(grep "^>" $TEMPPATHi_0 | wc -l)
  # append SAMPLEi to as many lines as there are sequences for the ith sample
  #[ $i = 1 ] && echo $(seq 1 $NUMSEQSi) | awk '{gsub("[0-9]+","'$SAMPLEi'\n")}1' | awk '{gsub(" ","")}1' > $SAMPLENAMES
  echo $(seq 1 $NUMSEQSi) | awk '{gsub("[0-9]+","'$SAMPLEi'\n")}1' | awk '{gsub(" ","")}1' | grep "\S" >> $SAMPLENAMES
  # append the new sequence names to NEWNAMES
  #[ $i = 1 ] && grep "^>" $TEMPPATHi | awk '{gsub("^>","")}1' > $NEWNAMES
  grep "^>" $TEMPPATHi | awk '{gsub("^>","")}1' >> $NEWNAMES
done

# paste the old and new sequence names and the associated sample names to the key file
paste $OLDNAMES $NEWNAMES $SAMPLENAMES >> $KEY

# save a three-column, space-delimited table (no header) with the new sample name, old sample name, and sequence type (sequence source) in the first, second, and third columns, respectively.
grep '_0000000001' $KEY | awk '{ print $2, $3}' | awk '{gsub("_0000000001","")}1' > $BLASTDBDIR"/allsamples_key_temp.txt"
echo $SEQTYPES | awk '{gsub(" ","\n")}1' > $BLASTDBDIR"/seqtypes_temp.txt"
paste -d " " $BLASTDBDIR"/allsamples_key_temp.txt" $BLASTDBDIR"/seqtypes_temp.txt" > $BLASTDBDIR"/allsamples_key.txt"
rm $BLASTDBDIR"/allsamples_key_temp.txt" $BLASTDBDIR"/seqtypes_temp.txt"

# Paths to the copied and updated files
TEMPPATHS=$(find $TEMPDIR -type f | sort)
## ##  # Run using pairwise approach or not
## ##  if [ "$PAIRWISE" = true ] ; then 

# directory where the actual blast databases should be saved
DATABASESDIR=$BLASTDBDIR'/databases'
[ ! -d "$DATABASESDIR" ] && mkdir $DATABASESDIR

# Prefix to use for each blast database
# DBPREFIX=$BLASTDBDIR'/blast_db' # <--- previous version
DBPREFIX=$DATABASESDIR'/blast_db'
# Make a blast database for each sample
for i in $(seq 1 $NUMSAMPLES); do
  PATHi=$(echo $TEMPPATHS | awk -v i=$i '{print $i}')
  DBSUFFIX=$(head -n 1 $PATHi | awk '{gsub("^>","")}1' | awk '{gsub("_.+","")}1')
  makeblastdb -in $PATHi -parse_seqids -dbtype nucl -out $DBPREFIX"_"$DBSUFFIX
done

## ##  else
## ##  	### concatenate all files listed in SEQFILE to 'allseqs.fa' in BLASTDBDIR
## ##  	ALLSEQS=$BLASTDBDIR'/allseqs.fa'
## ##  	[ ! -f $ALLSEQS ] && cat $TEMPPATHS > $ALLSEQS
## ##  	# Make the blast database
## ##  	makeblastdb -in $ALLSEQS -parse_seqids -dbtype nucl -out $BLASTDBDIR'/blast_db'
## ##  fi

#else
#	echo "Using existing database found in BLASTDBDIR"
#fi

###
# Run BLAST
###
# Parameters passed as arguments
## ## EVAL="${4:-'0.01'}"   # Expected value (e-value). Default is 0.01 (NCBI default is 10). Decreased EVAL to return fewer but better hits. Setting this too low might cause some homologous regions to be missed.
## ## MAXHSPs="${5:-'1'}"   # Maximum number of hits to return between a particular query-subject pair. NCBI default is 10, but here the default is 1.
## ## WORDSIZE="${6:-'11'}" # Minimum number of contiguous, identically matching base pairs to return a hit. Contiguous bases can be separated by gaps when dc-megablast is used; word_size must be 11 or 12 for dc-megablast, ...
## ## TASK="${7:-'blastn'}" # The algorithm to use. Default is 'blastn'. Other nucleotide options are 'dc-megablast', 'megablast', 'tblastn', and 'tblastx'.
# Parameters not passed as arguments
## ## GAPOPEN="5"         # Penalty multiplier for opening a gap. Increasing the value would likely return a set of hits with fewer gaps on avaerage, and possibly fewer total hits.
## ## GAPEXTEND="2"       # Penalty multiplier for increasing the length of a gap. Increasing the value would likely return a set of hits with shorter gaps on avaerage, and possibly fewer total hits.
## ## REWARD="2"          # reward value for each bp match
## ## PENALTY="-3"        # penalty value for each bp mismatch
## ## PERCIDENTITY="0"    # Percent identity cutoff. Minimum percent identity required between hits.
## ##  if [ "$PAIRWISE" = true ] ; then

# Directory that will hold blast hit tables
MATCHESDIR=$BLASTDBDIR'/matches'
[ ! -d "$MATCHESDIR" ] && mkdir $MATCHESDIR

for i in $(seq 1 $NUMSAMPLES); do
  for j in $(seq 1 $NUMSAMPLES); do
    # blast the sequences of the ith sample against the database of the jth sample, unless i=j
    if [ ! $i = $j ] ; then
      QUERYPATH=$(echo $TEMPPATHS | awk -v i=$i '{print $i}')
      QUERYi=$(head -n 1 $QUERYPATH | awk '{gsub("^>","")}1' | awk '{gsub("_.+","")}1')
      SUBJECTj=$(head -n 1 $(echo $TEMPPATHS | awk -v j=$j '{print $j}') | awk '{gsub("^>","")}1' | awk '{gsub("_.+","")}1')
      DBPATH=$DBPREFIX"_"$SUBJECTj
      # MATCHES=$BLASTDBDIR'/Q'$QUERYi'_S'$SUBJECTj'_matches.txt'
      MATCHES=$MATCHESDIR'/Q'$QUERYi'_S'$SUBJECTj'_matches.txt'
      # In the future blast_ij.sh will need to be in the same directory as 03.sh
      SHPATHij="/panfs/pfs.local/home/j926w878/work/conda/snakecap/blast_ij.sh"
      # sbatch --nodes=1 --ntasks-per-node=4 --mem=100Gb --time=24:00:00 --partition=bi $SHPATHij $TASK $DBPATH $QUERYPATH $MATCHES $MAXHSPs $EVAL 6 10 $WORDSIZE $GAPOPEN $GAPEXTEND $REWARD $PENALTY $PERCIDENTITY
      sbatch --nodes=1 --ntasks-per-node=4 --time=6:00:00 --partition=sixhour $SHPATHij $DBPATH $QUERYPATH $MATCHES
    fi
  done
done
## ##  else 
## ##  	MATCHES=$3 # Path where the output file (hit table in format 6) should be written.
## ##  	blastn -task $TASK -db $BLASTDBDIR'/blast_db' -query $ALLSEQS -out $MATCHES -max_hsps $MAXHSPs -evalue $EVAL -outfmt 6 -num_threads 10 -word_size $WORDSIZE -gapopen $GAPOPEN -gapextend $GAPEXTEND -reward $REWARD -penalty $PENALTY -perc_identity $PERCIDENTITY
## ##  done


