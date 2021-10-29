#!/bin/bash
# usage: sbatch --nodes=1 --mem=20Gb --partition=sixhour --time=6:00:00 /panfs/pfs.local/home/j926w878/work/conda/snakecap/test6_v2.sh "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10" "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/blast_database10/groupsMatrix_capVScap_capVStargets.txt"
#module load R
module load anaconda
conda activate py36

WORKDIR=$1
GROUPMATPATH=$2

# This produces the form that you want for the key file, except without the sequence names that dont belong to a group; need to get the ungrouped sequences by comparing to allseqs_key, and assign those to group0.
GROUPEDALIASPATH=$WORKDIR"/alias_groupedSeqs.txt"
cat $GROUPMATPATH | awk '{gsub("\t","\tgroup")}1' | awk '{gsub("\"","")}1' | awk 'NR>1' | awk '{print $1"\t"$1"_"$2}' | sort -V > $GROUPEDALIASPATH

### Compares newnames.txt to $GROUPEDALIASPATH{$1} in a way like setdiff() in R, but uses the comm bash command
# creates a file with values from the first column of $GROUPMATPATH, removes double quotes, and removes the header line; the sort -V option behaves like the gtools::mixedsort R function, but shouldnt be used until after using comm
# cat $GROUPMATPATH | awk '{gsub("\"","")}1' | awk 'NR>1' | awk '{print $1}' | sort -V > $WORKDIR"/groupedSeqs.txt"
cat $GROUPMATPATH | awk '{gsub("\"","")}1' | awk 'NR>1' | awk '{print $1}' | sort > $WORKDIR"/groupedSeqs.txt"
# creates a file with the set of unique samples in $GROUPMATPATH; keeps the underscore so that partial matches (e.g., sample1 and sample10) do not cause confusion.
cat $GROUPMATPATH | awk '{gsub("\"","")}1' | awk 'NR>1' | awk '{print $1}' | awk '{gsub("_.+","_")}1' | sort -V | uniq > $WORKDIR"/samplesWithGroupedSeqs.txt"
# creates a file with the lines of newnames.txt for samples present in GROUPMATPATH.
grep -f $WORKDIR"/samplesWithGroupedSeqs.txt" $WORKDIR"/newnames.txt" | sort > $WORKDIR"/AllSequenceNamesForSamplesWithGroupedSeqs.txt"
# creates a file with the names of the ungrouped sequences on seperate lines for the samples with at least one grouped sequence; -3 flag means to suppress lines shared in both of the input files; --check-order means to check that the lines are correctly sorted, --nocheck-order means do not check the order, which is what we want because we sorted in different way than expected by comm
comm -3 --check-order $WORKDIR"/AllSequenceNamesForSamplesWithGroupedSeqs.txt" $WORKDIR"/groupedSeqs.txt" | sort -V > $WORKDIR"/ungroupedSeqs.txt"
### The sum of the number of lines in ungroupedSeqs and groupedSeqs is equal the number of lines in AllSequenceNamesForSamplesWithGroupedSeqs.txt, indicating that comm did what we wanted

# creates a file with a structure like $ALIASPATH, but only includes the ungrouped sequences; "_ungrouped" will be appended to the names.
UNGROUPEDALIASPATH=$WORKDIR"/alias_ungroupedSeqs.txt"
cat $WORKDIR"/ungroupedSeqs.txt" | awk '{print $1"\t"$1"_ungrouped"}' > $UNGROUPEDALIASPATH
# concatenate $ALIASPATH and $UNGROUPEDALIASPATH
cat $GROUPEDALIASPATH $UNGROUPEDALIASPATH | sort -V > $WORKDIR"/alias.txt"

# Use seqkit to create a copy of the set of fasta files in tempseqs for the samples with at least one grouped sequence, but replace sequence IDs using the names in column two of alias.txt
GROUPSEQSDIR=$WORKDIR"/groupseqs"
[ ! -d "GROUPSEQSDIR" ] && mkdir $GROUPSEQSDIR
TEMPDIR="$WORKDIR""/tempseqs"
SAMPLESKEY=$WORKDIR"/allsamples_key.txt"
SAMPLES=$(cat $WORKDIR"/samplesWithGroupedSeqs.txt" | awk '{gsub("_","")}1')
NUMSAMPLES=$(echo $SAMPLES | wc -w)

for VARi in $(seq 1 $NUMSAMPLES) ; do
  SAMPLEi=$(echo $SAMPLES | awk -v i="$VARi" '{print $i}')
  NAMEi=$(grep -w $SAMPLEi $SAMPLESKEY | awk '{print $2}')
  SEQPATHi=$(find $TEMPDIR -name "*"$NAMEi"*")
  # define path to a temporary fasta file that will hold the fasta with replaced sequence IDs for the ith sample
  TEMPFASTA=$WORKDIR"/groupseqs/"$SAMPLEi".fa"
  # Now use seqkit
  seqkit replace -p '(.+)$' -r '{kv}' -k $WORKDIR"/alias.txt" $SEQPATHi -w 0 > $TEMPFASTA
  #seqkit split $TEMPFASTA -i --id-regexp "^sample[0-9]+_[0-9]+_([\w]+[0-9]+)$" -2
  seqkit split $TEMPFASTA -i --id-regexp "^.+[0-9]_[0-9]+_([\w]+)$" -2 -w 0
  # remove TEMPFASTA
  rm $TEMPFASTA
  # Im not sure if I should remove the .fai file, so wont for now
done

# Groups names (1, 2, ..., N; ungrouped)
GROUPNAMES=$(awk -F '\t' 'NR>1{print "group"$2}' $GROUPMATPATH | sort | uniq | sort -V)" ungrouped"
# For each group, gather all sequences of the group and put them in a fasta file.
NGROUPS=$(echo $GROUPNAMES | wc -w)
for VARj in $(seq 1 $NGROUPS) ; do
  GROUPNAMEi=$(echo $GROUPNAMES | awk -v i=$VARj '{print $i}')
  cat $(find $GROUPSEQSDIR -regex '.*[^/]*_'${GROUPNAMEi}'.fasta') > $GROUPSEQSDIR"/"$GROUPNAMEi".fasta"
done

