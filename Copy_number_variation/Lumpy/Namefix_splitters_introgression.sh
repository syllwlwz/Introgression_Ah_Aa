#!/bin/bash

#$ -t 1-178
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Merge
#$ -l vf=56G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat /prj/pflaphy-gscan/Final_samples_introgression_single.list | head -n $SGE_TASK_ID | tail -n 1)
java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/$file.splitters O=/prj/pflaphy-gscan/dedup/$file.namefixed.splitters SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=$file RGPL=illumina RGPU=AAAAAA RGSM=$file\ 
 2> /prj/pflaphy-gscan/dedup/$file.namefixed.splitters.err
