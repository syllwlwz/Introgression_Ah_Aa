#!/bin/bash

#Alignment

#$ -t 1-7
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Alignment_lyrata
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 16

file=$(head -n $SGE_TASK_ID Introgression_samples_mapping3.list | tail -n 1)
dir=$(head -n $SGE_TASK_ID Introgression_samples_mapping3.list | tail -n 1)
PU=$(head -n $SGE_TASK_ID Introgression_samples_mapping3.list | tail -n 1)

./software/bwa-0.7.17/bwa mem -t 16 -R $(echo "@RG\tID:$file\tSM:$dir\tLB:$dir\tPU:$PU\tPL:ILLUMINA") -k 10 Alyrata_107.fa -M cutadapt_introgression/$file.pass.1.fastq.sub.gz cutadapt_introgression/$file.pass.2.fastq.sub.gz > aligned/$file.sam 2> aligned/$file.sam.err

#t: number of threads
#k: minimum length of seed region
 
./software/samtools-1.11/samtools view -b aligned/$file.sam | ./software/samtools-1.11/samtools sort -T $file > aligned/$file.sort.bam 2> aligned/$file.sort.bam.err

#b: output in bam format
#T: write temporary files

./software/samtools-1.11/samtools index aligned/$file.sort.bam

./software/samtools-1.11/samtools idxstats aligned/$file.sort.bam > aligned/$file.idxstats

./software/samtools-1.11/samtools flagstat aligned/$file.sort.bam > aligned/$file.flagstat

rm aligned/$file.sam

