#!/bin/bash

#cutadapt version 1.11

#$ -t 31-38
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Fastqc_outgroup
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 6

file=$(cat Outgroup.list | head -n $SGE_TASK_ID | tail -n 1)
./software/FastQC/fastqc -t 6 $file.sub.fastq.gz -o Thaliana_1001genomes/ 2> Thaliana_1001genomes/Fastqc_$SGE_TASK_ID.err

