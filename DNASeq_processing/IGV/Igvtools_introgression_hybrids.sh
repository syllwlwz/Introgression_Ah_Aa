##!/bin/bash

#$ -t 1-11
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N igvtools_introgression
#$ -l vf=8G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Concat_alignment.list | head -n $SGE_TASK_ID | tail -n 1)

./software/IGV_2.11.2/igvtools count --minMapQuality 0 -w 100 --includeDuplicates /prj/pflaphy-gscan/dedup/$file.dedup.bam /prj/pflaphy-gscan/dedup/$file.dedup.bam.MQ0.tdf /prj/pflaphy-gscan/Halleri_arenosa.genome\
 2> ./software/IGV_2.11.2/Tdf.MQ0.$file.err
