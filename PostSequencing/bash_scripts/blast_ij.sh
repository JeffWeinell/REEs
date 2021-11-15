#!/bin/bash
module load 'blast+/2.9.0'
#blastn -task $1 -db $2 -query $3 -out $4 -max_hsps $5 -evalue $6 -outfmt $7 -num_threads $8 -word_size $9 -gapopen $10 -gapextend $11 -reward $12 -penalty $13 -perc_identity $14
#blastn -task 'blastn' -db "$1" -query "$2" -out "$3" -max_hsps 1 -evalue 0.01 -outfmt 6 -num_threads 10 -word_size 10 -gapopen 5 -gapextend 2 -reward 2 -penalty '-3' -perc_identity 0
echo $1"_"$2
blastn -task 'blastn' -db "$1" -query "$2" -out "$3" -max_hsps 1 -evalue 0.0000000001 -outfmt 6 -num_threads 10 -word_size 100 -gapopen 5 -gapextend 2 -reward 2 -penalty '-3' -perc_identity 0
