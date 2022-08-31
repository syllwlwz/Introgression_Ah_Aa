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

java -Xmx50g -jar /prj/pflaphy-gscan/software/picard.jar MarkDuplicates I=/prj/pflaphy-gscan/aligned/Kowa_003_h02_GAGATTCC-TATAGCCT.sort.bam I=/prj/pflaphy-gscan/aligned/Kowa_h02.sort.bam O=/prj/pflaphy-gscan/dedup/Kowa_h02comb.dedup.bam M=/prj/pflaphy-gscan/dedup/duplicateMetricsFile_Kowa_h02comb VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true 2> /prj/pflaphy-gscan/dedup/Kowa_h02comb.err
