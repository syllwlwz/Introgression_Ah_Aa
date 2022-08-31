#!/bin/bash

#run from pflaphy-gscan

#$ -t 1-187
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N discordant_RP
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Final_samples_introgression.list | head -n $SGE_TASK_ID | tail -n 1)
# Extract the discordant paired-end alignments.
./software/samtools-1.11/samtools view -b -F 1294 dedup/$file.dedup.bam > dedup/$file.discordants.unsorted.bam

# Sort alignments
./software/samtools-1.11/samtools sort -T $file dedup/$file.discordants.unsorted.bam > dedup/$file.discordants
