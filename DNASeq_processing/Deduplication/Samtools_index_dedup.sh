#!/bin/bash

#run from pflaphy-gscan

#$ -t 201-204
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N samtools_index
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Total_samples_introgression.list | head -n $SGE_TASK_ID | tail -n 1)
./software/samtools-1.11/samtools index dedup/$file.dedup.bam 2> dedup/$file.index.err

