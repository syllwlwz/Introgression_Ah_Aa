#!/bin/bash

#cutadapt version 2.10

#$ -t 1-25
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cutadapt
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

file=$(cat Samples_Vinod/Samples_Vinod.list | head -n $SGE_TASK_ID | tail -n 1)
./software/cutadapt -e 0.15 -O 4 -m 120 --nextseq-trim=20 -j 12 -q 20 --max-n 20 -a "G{150}" -A "G{150}" -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA' -n 2 -o cutadapt_introgression/$file.R1.cutadapt.fastq.gz -p cutadapt_introgression/$file.R2.cutadapt.fastq.gz Samples_Vinod/$file.R1.fq.gz Samples_Vinod/$file.R2.fq.gz > cutadapt_introgression/$file.cutadapt 2> cutadapt_introgression/$file.cutadapt.err

#chmod 755 cutadapt.sh
#nohup ./cutadapt.sh &


