#!/bin/bash

#Deduplication
#picard tools 2.21.9

#$ -t 1-25
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Deduplication
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(head -n $SGE_TASK_ID Samples_Vinod/Samples_Vinod.list | tail -n 1)

java -Xmx50g -jar /prj/pflaphy-gscan/software/picard.jar MarkDuplicates I=/prj/pflaphy-gscan/aligned/$file.sort.bam O=/prj/pflaphy-gscan/dedup/$file.dedup.bam M=/prj/pflaphy-gscan/dedup/duplicateMetricsFile_$file VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=false OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 REMOVE_SEQUENCING_DUPLICATES=true 2> /prj/pflaphy-gscan/dedup/$file.err
