#!/bin/bash

#$ -t 1-60
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Bam_sort
#$ -l vf=7G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Cutadapt_introrna.list | head -n $SGE_TASK_ID | tail -n 1)
software/samtools-1.14/samtools sort Error_corrected2/$file'_outfile.bam' -o Error_corrected2/$file.sorted.bam

