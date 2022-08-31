##!/bin/bash

#$ -t 1-187
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N igvtools_introgression
#$ -l vf=8G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Igvtools.list | head -n $SGE_TASK_ID | tail -n 1)
#HMA4
#./software/IGV_2.11.2/igvtools count -w 1 --minMapQuality 25 --includeDuplicates /prj/pflaphy-gscan/dedup/$file /prj/pflaphy-gscan/dedup/$file.tdf /prj/pflaphy-gscan/A_lyrata.genome\
# 2> ./software/IGV_2.11.2/Tdf.$file.err

./software/IGV_2.11.2/igvtools count --minMapQuality 25 -w 100 --includeDuplicates /prj/pflaphy-gscan/dedup/$file /prj/pflaphy-gscan/dedup/$file.tdf /prj/pflaphy-gscan/A_lyrata.genome\
 2> ./software/IGV_2.11.2/Tdf.$file.err
