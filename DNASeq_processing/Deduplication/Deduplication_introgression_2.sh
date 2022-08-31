#!/bin/bash

#Deduplication
#picard tools 2.21.9

#$ -t 1-27
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Deduplication
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(head -n $SGE_TASK_ID Outgroup_dedup.list | tail -n 1)
out=$(echo $file | sed 's/\.1//g')

java -Xmx50g -jar /prj/pflaphy-gscan/software/picard.jar MarkDuplicates I=/prj/pflaphy-gscan/aligned/$file.sort.bam O=/prj/pflaphy-gscan/dedup/$out.dedup.bam M=/prj/pflaphy-gscan/dedup/duplicateMetricsFile_$file VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true 2> /prj/pflaphy-gscan/dedup/$file.err
