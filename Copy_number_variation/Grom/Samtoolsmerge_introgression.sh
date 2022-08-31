#!/bin/bash

#samtools 1.11.1 with htslib 1.11.1

#$ -t 1-6
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N samtools_merge
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Grom_introgression.list | head -n $SGE_TASK_ID | tail -n 1)
./software/samtools-1.11/samtools merge -b $file'_halleri_full.list' GROM_analysis/$file'_halleri.bam' 2> GROM_analysis/$file'_halleri.err'

./software/samtools-1.11/samtools merge -b $file'_arenosa_full.list' GROM_analysis/$file'_arenosa.bam' 2> GROM_analysis/$file'_arenosa.err'

