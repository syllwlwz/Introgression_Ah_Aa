#!/bin/bash

#cutadapt version 2.10

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cutadapt
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

./software/cutadapt -e 0.15 -O 4 -m 50 -j 12 -q 20 --max-n 20 -A "G{150}" -a "G{150}" -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA' -n 2 -o Thaliana_1001genomes/SRR7637309.R1.cutadapt.fastq.gz -p Thaliana_1001genomes/SRR7637309.R2.cutadapt.fastq.gz Thaliana_1001genomes/SRR7637309_1.fastq.gz Thaliana_1001genomes/SRR7637309_2.fastq.gz > Thaliana_1001genomes/SRR7637309.cutadapt 2> Thaliana_1001genomes/SRR7637309.cutadapt.err

#chmod 755 cutadapt.sh
#nohup ./cutadapt.sh &


