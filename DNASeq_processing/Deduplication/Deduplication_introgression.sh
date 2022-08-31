#!/bin/bash

#Deduplication
#picard tools 2.21.9

#$ -t 1-96
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Deduplication
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(head -n $SGE_TASK_ID filelist_introgression.list | tail -n 1)
files=$(grep $file'_' Introgression_raw_data.list | cut -f 4 -d "/" | sort | uniq | sed -e 's/$/.sort.bam/' | sed -e 's/^/I=\/prj\/pflaphy-gscan\/aligned\//' | tr "\n" "\t")
echo $file
echo $files

java -Xmx50g -jar /prj/pflaphy-gscan/software/picard.jar MarkDuplicates $files O=/prj/pflaphy-gscan/dedup/$file.dedup.bam M=/prj/pflaphy-gscan/dedup/duplicateMetricsFile_$file VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=false OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 REMOVE_SEQUENCING_DUPLICATES=true 2> /prj/pflaphy-gscan/dedup/$file.err
