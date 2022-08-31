#!/bin/bash

#run from pflaphy-gscan

#$ -t 1-187
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Distrib_lumpy
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Final_samples_introgression.list | head -n $SGE_TASK_ID | tail -n 1)
RL=$(software/samtools-1.11/samtools view dedup/$file.dedup.bam | awk '{print length($10)}' | head -1000 | sort -nu | tail -n 1)

./software/samtools-1.11/samtools view dedup/$file.dedup.bam | tail -n+100000 | software/lumpy-sv/scripts/pairend_distro.py -r $RL -X 4 -N 10000 -o dedup/$file.lib1.histo > dedup/$file.distrib
