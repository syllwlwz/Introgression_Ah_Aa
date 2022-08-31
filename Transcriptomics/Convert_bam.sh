#!/bin/bash

#$ -t 1-60
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Bam_convert
#$ -l vf=7G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Cutadapt_introrna.list | head -n $SGE_TASK_ID | tail -n 1)
software/samtools-1.14/samtools view -bo mapped_lyrata2/$file.sorted.bam mapped_lyrata2/$file'.sam' 2> mapped_lyrata2/$file.convert.err

