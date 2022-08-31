#!/bin/bash

#GROM
#samtools 1.3.1 with htslib 1.3.1
#run from pflaphy-gscan

#$ -t 1-6
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N samtools_index
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Grom_introgression.list | head -n $SGE_TASK_ID | tail -n 1)
./software/samtools-1.11/samtools index GROM_analysis/$file'_halleri.sorted.bam' GROM_analysis/$file'_halleri.sorted.bam.bai' 2> GROM_analysis/$file'_halleri.err'
./software/samtools-1.11/samtools index GROM_analysis/$file'_arenosa.sorted.bam' GROM_analysis/$file'_arenosa.sorted.bam.bai' 2> GROM_analysis/$file'_arenosa.err'

