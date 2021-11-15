#!/bin/bash
# usage: This is a subprogram that is called by test6_v2.sh
module load anaconda
conda activate py36

# variables passed from test6_v2.sh
WORKDIR=$1
GROUPSEQSDIR=$2
SAMPLEi=$3
SAMPLESKEY=$4
SEQPATHi=$5

# Name to use for the temporary fasta file
TEMPFASTA=$GROUPSEQSDIR"/"$SAMPLEi".fa"
# Use seqkit to append the group IDs to sequence IDs
seqkit replace -p '(.+)$' -r '{kv}' -k $WORKDIR"/alias.txt" $SEQPATHi -w 0 > $TEMPFASTA
# For a particular sample, split sequences into separate fasta files by group
seqkit split $TEMPFASTA -i --id-regexp "^.+[0-9]_[0-9]+_([\w]+)$" -2 -w 0
# remove TEMPFASTA
rm $TEMPFASTA
