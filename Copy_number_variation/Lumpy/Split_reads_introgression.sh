#!/bin/bash

#run from pflaphy-gscan

#$ -t 1-187
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Split_reads
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Final_samples_introgression.list | head -n $SGE_TASK_ID | tail -n 1)

# Extract the split-read alignments
./software/samtools-1.11/samtools view -h dedup/$file.dedup.bam | software/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | ./software/samtools-1.11/samtools view -Sb - > dedup/$file.splitters.unsorted.bam

# Sort alignments
./software/samtools-1.11/samtools sort -T $file dedup/$file.splitters.unsorted.bam > dedup/$file.splitters
