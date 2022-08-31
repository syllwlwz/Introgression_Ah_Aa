#!/bin/bash

#run from pflaphy-gscan

#$ -t 1-178
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N samtools_index
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat /prj/pflaphy-gscan/Final_samples_introgression_single.list | head -n $SGE_TASK_ID | tail -n 1)
./software/samtools-1.11/samtools index dedup/$file.namefixed.dedup.bam dedup/$file.namefixed.dedup.bam.bai 2> dedup/$file.namefixed.dedup.index.err

