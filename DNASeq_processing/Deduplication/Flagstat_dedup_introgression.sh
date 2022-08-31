#!/bin/bash

#Alignment

#$ -t 201-204
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Sum_dedup_lyrata
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(head -n $SGE_TASK_ID Total_samples_introgression.list | tail -n 1)

./software/samtools-1.11/samtools idxstats dedup/$file.dedup.bam > dedup/$file.idxstats

./software/samtools-1.11/samtools flagstat dedup/$file.dedup.bam > dedup/$file.flagstat

