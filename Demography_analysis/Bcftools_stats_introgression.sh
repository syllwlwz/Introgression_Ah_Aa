#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Bcftools_stats
#$ -l vf=200M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

./software/bcftools-1.12/bcftools stats -S Introgression_samples_final_vcf.list Structure/LD_Pruned/All_introgression_4dgsites.LD_Pruned.vcf --threads 4 > Structure/LD_Pruned/All_introgression_4dgsites.LD_Pruned.stats.vcf 2> Structure/LD_Pruned/All_introgression_4dgsites.LD_Pruned.stats.err

