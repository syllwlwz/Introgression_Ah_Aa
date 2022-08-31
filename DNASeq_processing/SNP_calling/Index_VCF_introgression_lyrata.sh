#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N IndexVCF
#$ -l vf=56G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar IndexFeatureFile -I /prj/pflaphy-gscan/HC2_introgression/Introgression_lyrata.full.vcf.gz  2> /prj/pflaphy-gscan/HC2_introgression/Lyrata_Indexcvcf_introgression.err 
