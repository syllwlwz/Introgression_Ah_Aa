##!/bin/bash

#$ -t 44
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N igvtools_introgression
#$ -l vf=8G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(head -n $SGE_TASK_ID Dedup_BAC.list | tail -n 1)
#file=$(head -n $SGE_TASK_ID aligned_BACs/Comb.list | tail -n 1)
./software/IGV_2.11.2/igvtools count --minMapQuality 25 -w 1 --includeDuplicates /prj/pflaphy-gscan/aligned_BACs/$file.dedup.bam /prj/pflaphy-gscan/aligned_BACs/$file.tdf /prj/pflaphy-gscan/HMA4_BACs.genome\
 2> ./software/IGV_2.11.2/BACs.Tdf.$file.err
