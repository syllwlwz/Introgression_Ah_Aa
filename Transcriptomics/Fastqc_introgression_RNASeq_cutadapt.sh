#!/bin/bash

#$ -t 1-60
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Fastqc_introgression2
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 6

file=$(cat Sample.list | head -n $SGE_TASK_ID | tail -n 1)
./software/FastQC/fastqc -t 6 cutadapt/$file.1.cutadapt.fastq.gz -o Fastqc_results_introgression_cutadapt/ 2> Fastqc_results_introgression_cutadapt/Fastqc_cutadapt1_$SGE_TASK_ID.err
./software/FastQC/fastqc -t 6 cutadapt/$file.2.cutadapt.fastq.gz -o Fastqc_results_introgression_cutadapt/ 2> Fastqc_results_introgression_cutadapt/Fastqc_cutadapt2_$SGE_TASK_ID.err

