#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N GatherVCF
#$ -l vf=40G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -jar software/picard.jar GatherVcfs I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_1.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_2.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_3.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_4.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_5.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_6.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_7.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_8.vcf.gz I=/prj/pflaphy-gscan/HC2_introgression/HC2_introgression_scaffold_9.vcf.gz O=/prj/pflaphy-gscan/HC2_introgression/Introgression.vcf.gz  2> /prj/pflaphy-gscan/HC2_introgression/Gathervcfs_introgression.err 
