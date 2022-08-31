#!/bin/bash

#GATK4.1.9.0

#$ -t 1-9
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N GVCF
#$ -l vf=24G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

Scaffold=$(head -n $SGE_TASK_ID Mainscaffolds.list | tail -n 1)

java -Xmx80g -Xms80g -XX:ConcGCThreads=4 -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenotypeGVCFs -R /prj/pflaphy-gscan/Alyrata_107.fa -V gendb:///prj/pflaphy-gscan/DB_introgression_lyrata_full_$Scaffold/ -L $Scaffold --heterozygosity 0.005\
 --indel-heterozygosity 0.001 -O /prj/pflaphy-gscan/HC2_introgression/HC2_introgression_lyrata_full_$Scaffold.vcf.gz --tmp-dir /prj/pflaphy-gscan/HC2_temp/ -DF NotDuplicateReadFilter --genomicsdb-shared-posixfs-optimizations true\
 --include-non-variant-sites 2> /prj/pflaphy-gscan/HC2_introgression/HC2_introgression_lyrata_full_$Scaffold.err 
