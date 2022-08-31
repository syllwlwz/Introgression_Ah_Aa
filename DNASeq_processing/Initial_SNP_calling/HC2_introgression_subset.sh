#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N HC2_introgression
#$ -l vf=15G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenotypeGVCFs -R /prj/pflaphy-gscan/Alyrata_107.fa -V /prj/pflaphy-gscan/HC2_introgression/Introgression_subset.g.vcf.gz -O /prj/pflaphy-gscan/HC2_introgression/HC2_introgression_subset.vcf.gz --indel-heterozygosity 0.01 --heterozygosity 0.02 2> /prj/pflaphy-gscan/HC2_introgression/HC2_introgression_subset.err 

