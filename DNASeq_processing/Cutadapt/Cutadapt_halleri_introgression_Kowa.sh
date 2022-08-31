#!/bin/bash

#$ -t 1-7
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cutadapt_halleri
#$ -l vf=600M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

file=$(cat Kowa.list | head -n $SGE_TASK_ID | tail -n 1)
./software/cutadapt -e 0.15 -O 4 -m 50 -j 12 --max-n 20 -n 2 -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA' -o cutadapt_introgression/$file.R1.cutadapt.fastq.gz -p cutadapt_introgression/$file.R2.cutadapt.fastq.gz Fastq_files/$file.R1.fastq.gz Fastq_files/$file.R2.fastq.gz > cutadapt_introgression/$file.cutadapt 2> cutadapt_introgression/$file.cutadapt.err

