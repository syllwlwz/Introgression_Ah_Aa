#!/bin/bash

#GROM

#$ -t 1-6
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N GROM
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

file=$(cat Grom_introgression.list | head -n $SGE_TASK_ID | tail -n 1)
./software/GROM/dist/GROM -P 6 -b 25 -q 25 -v 0.01 -e 0.01 -V 0.01 -U 4 -A 4 -p 2 -s 5 -i GROM_analysis/$file'_halleri.sorted.bam' -r Alyrata_107.fa -o GROM_analysis/$file.halleri.predictions.vcf > GROM_analysis/$file.halleri.out 2> GROM_analysis/$file.halleri.GROM.err

#Select INDELs and CNVs
grep "#" GROM_analysis/$file.halleri.predictions.vcf > GROM_analysis/$file.halleri.cnvs.vcf
grep "SD:Z:CN:CS"  GROM_analysis/$file.halleri.predictions.vcf >> GROM_analysis/$file.halleri.cnvs.vcf

