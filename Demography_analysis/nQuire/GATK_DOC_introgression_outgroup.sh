#!/bin/bash

#GATK version 3.7
#cohorts separately

#$ -t 26
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Average_coverage_introgression
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

file=$(cat Outgroup_renamed.list | head -n $SGE_TASK_ID | tail -n 1)

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar DepthOfCoverage -I dedup/$file.dedup.bam -R Alyrata_107.fa -O dedup/$file.lyrata.mean_DOC_all --min-base-quality  25 --omit-depth-output-at-each-base --omit-interval-statistics --omit-locus-table -DF NotDuplicateReadFilter -L Mainscaffolds2.bed 2> dedup/$file.lyrata.mean_DOC.err

