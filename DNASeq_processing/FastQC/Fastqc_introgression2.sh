#!/bin/bash

#cutadapt version 1.11

#$ -t 1-192
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Fastqc_introgression2
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 6

file=$(cat Fastqc_introgression2.list | head -n $SGE_TASK_ID | tail -n 1)
./software/FastQC/fastqc -t 6 $file -o Fastqc_results_introgression2/ 2> Fastqc_results_introgression2/Fastqc2_$SGE_TASK_ID.err

