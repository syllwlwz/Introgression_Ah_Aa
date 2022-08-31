#!/bin/bash

#cutadapt version 3.5

#$ -t 27
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cutadapt_rna
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Cutadapt_introrna.list | head -n $SGE_TASK_ID | tail -n 1)
./software/cutadapt -e 0.15 -O 4 -m 120 --nextseq-trim=20 -q 20 -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' -a "A{100}" -A "A{100}" -a "G{100}" -A "G{100}" -n 2 -o cutadapt/$file.1.cutadapt.fastq.gz -p cutadapt/$file.2.cutadapt.fastq.gz\
 Vero_RNASeq_intro/X204SC21090574-Z01-F001_02/raw_data/$file/$file'_1.fq.gz' Vero_RNASeq_intro/X204SC21090574-Z01-F001_02/raw_data/$file/$file'_2.fq.gz' > cutadapt/$file.cutadapt 2> cutadapt/$file.cutadapt.err

