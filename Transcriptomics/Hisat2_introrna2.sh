#!/bin/bash

#$ -t 1-60
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Hisat2_lyrata
#$ -l vf=200M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

file=$(cat Cutadapt_introrna.list | head -n $SGE_TASK_ID | tail -n 1)
./software/hisat2-2.2.0/hisat2 -p 12 --no-unal -q -x Alyrata_all.index -1 cutadapt/$file.1.cutadapt.fastq.gz -2 cutadapt/$file.2.cutadapt.fastq.gz -S mapped_lyrata2/$file'.sam' --summary-file mapped_lyrata2/$file'.summary.txt' --rna-strandness RF --max-intronlen 20000 -L 10 --pen-noncansplice 16 --mp 2,0 --score-min L,0,-0.4 --rdg 3,1 --rfg 3,1 --pen-canintronlen G,-0.5,0.1 -R 3 -i S,1,0.50 --bowtie2-dp 2 -k 50 2> mapped_lyrata2/$file'.err'
