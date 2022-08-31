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

java -Xmx50g -jar /prj/pflaphy-gscan/software/picard.jar MarkDuplicates I=/prj/pflaphy-gscan/aligned/Zako001h01.sort.bam I=/prj/pflaphy-gscan/aligned/Zako_h09.sort.bam O=/prj/pflaphy-gscan/dedup/Zako_h09comb.dedup.bam M=/prj/pflaphy-gscan/dedup/duplicateMetricsFile_Zako_h09comb VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 2> /prj/pflaphy-gscan/dedup/Zako_h09comb.err
