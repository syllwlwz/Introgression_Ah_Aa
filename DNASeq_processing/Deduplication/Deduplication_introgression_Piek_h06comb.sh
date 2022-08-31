#!/bin/bash

#Deduplication
#picard tools 2.21.9

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Deduplication
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

files=$(grep L37'_' Introgression_raw_data.list | cut -f 4 -d "/" | sort | uniq | sed -e 's/$/.sort.bam/' | sed -e 's/^/I=\/prj\/pflaphy-gscan\/aligned\//' | tr "\n" "\t")
echo $files

java -Xmx50g -jar /prj/pflaphy-gscan/software/picard.jar MarkDuplicates $files I=/prj/pflaphy-gscan/aligned/Piek6.sort.bam O=/prj/pflaphy-gscan/dedup/Piek_h06comb.dedup.bam M=/prj/pflaphy-gscan/dedup/duplicateMetricsFile_Piek_h06comb VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 2> /prj/pflaphy-gscan/dedup/Piek_h06comb.err
