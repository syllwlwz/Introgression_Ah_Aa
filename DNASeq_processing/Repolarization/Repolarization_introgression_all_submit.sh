#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Repol_all2
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#./software/R-4.1.0/bin/Rscript Repolarization_introgression_all.R
#sed -i 's/NA/./g' Filtered_introgression2/All2_repol.vcf
grep "#" /prj/pflaphy-gscan/Filtered_introgression2/All2.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/All2_repol.vcf > /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf

