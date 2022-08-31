#!/bin/bash

#run from pflaphy-gscan

#$ -t 1-9
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N DBimport
#$ -l vf=25G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

Scaffold=$(head -n $SGE_TASK_ID Mainscaffolds.list | tail -n 1)

java -Xmx80g -Xms80g -XX:ConcGCThreads=4 -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport -R /prj/pflaphy-gscan/Alyrata_107.fa\
 -V /prj/pflaphy-gscan/HC1_introgression/A_ce_1.full.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/A_ce_2.full.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/A_ce_3.full.g.vcf.gz\
 --use-jdk-inflater --tmp-dir /prj/pflaphy-gscan/HC2_temp/ -DF NotDuplicateReadFilter -L $Scaffold --genomicsdb-workspace-path /prj/pflaphy-gscan/DB_introgression_cebennensis_$Scaffold/ --batch-size 50 --max-num-intervals-to-import-in-parallel 4 --merge-input-intervals true --genomicsdb-shared-posixfs-optimizations true 2> /prj/pflaphy-gscan/HC2_introgression/GenomicsDBimport_cebennensis_$Scaffold.err


