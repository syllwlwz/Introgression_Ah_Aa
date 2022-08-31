#!/bin/bash

#cutadapt version 2.10

#$ -t 1-768
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cutadapt
#$ -l vf=600M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

file=$(cat Introgression_raw_data.list | head -n $SGE_TASK_ID | tail -n 1)
dir=$(head -n $SGE_TASK_ID Introgression_raw_data.list | tail -n 1 | cut -f 4 -d "/")
./software/cutadapt -e 0.15 -O 4 -m 120 --nextseq-trim=20 -j 12 -q 20 --max-n 20 -a "G{150}" -A "G{150}" -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA' -n 2 -o cutadapt_introgression/$dir.R1.cutadapt.fastq.gz -p cutadapt_introgression/$dir.R2.cutadapt.fastq.gz $file'_1.fq.gz' $file'_1.fq.gz' > cutadapt_introgression/$dir.cutadapt 2> cutadapt_introgression/$dir.cutadapt.err

#chmod 755 cutadapt.sh
#nohup ./cutadapt.sh &


