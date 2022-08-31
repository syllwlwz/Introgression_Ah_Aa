#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Vcftogenotype
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

python3 software/genomics_general/VCF_processing/parseVCF.py -i  /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr_fil.vcf.gz --skipIndels --skipMono --ploidyFile Twisst/MiashalleriKowaarenosa_Vcftogeno_groups.txt | gzip > Twisst/MiashalleriKowaarenosalyrata_SNPs_fil2.geno.gz 2> Twisst/MiashalleriKowaarenosalyrata_SNPs_fil2.geno.err
