#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Bcftools_rename
#$ -l vf=200M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

./software/bcftools-1.12/bcftools reheader -s Newvcfnames_halleri.txt HC2_introgression/Introgression_halleri.vcf.gz -o HC2_introgression/Introgression_halleri_renamed.vcf.gz --threads 4 2> HC2_introgression/Introgression_halleri_renamed.err

./software/bcftools-1.12/bcftools query -l HC2_introgression/Introgression_halleri_renamed.vcf.gz | sort > Newvcforder.txt
./software/bcftools-1.12/bcftools view -S Newvcforder.txt HC2_introgression/Introgression_halleri_renamed.vcf.gz -o HC2_introgression/Introgression_halleri_renamed_sorted.vcf.gz --threads 4 2> HC2_introgression/Introgression_halleri_renamed_sorted.err

