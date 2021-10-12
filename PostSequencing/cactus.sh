#!/bin/bash
# usage: cactus.sh </path/to/input/seqFile> </path/to/output.HAL> [/path/to/jobstore] [/path/to/working/directory]
# sbatch --nodes=1 --ntasks-per-node=10 --mem=200Gb --time=120:00:00 --partition=bi '/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/cactus.sh' '/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/seqFile.txt' '/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Cactus/alignment.HAL' 'jobStore'
module use /home/j926w878/sw/modules
module load cactus

SEQFILE=$1
HAL=$2
JOBSTORE=${3:-"jobStore"}
WORKDIR=${4:-"PWD"} # WORKDIR must complete the paths the sequence files in SEQFILE

# SEQPATHS=$(awk -F "\"* \"*" '!/^($|#)/NR>1{print $2}' $SEQFILE)

# Go here for help: 'https://github.com/ComparativeGenomicsToolkit/cactus'
# Cactus usage in general: cactus <jobStore> <seqFile> <outputHAL> [options]

# cd "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples"
cd $WORKDIR
cactus $JOBSTORE $SEQFILE $HAL --root sr --binariesMode local
