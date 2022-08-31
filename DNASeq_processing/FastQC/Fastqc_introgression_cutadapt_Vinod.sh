#!/bin/bash

#cutadapt version 1.11

#$ -t 1-50
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Fastqc_introgression
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 6

file=$(cat Samples_Vinod/Samples_Vinod_total.list | head -n $SGE_TASK_ID | tail -n 1)
./software/FastQC/fastqc -t 6 cutadapt_introgression/$file -o Fastqc_results_introgression/ 2> Fastqc_results_introgression/Fastqc_$SGE_TASK_ID.err

