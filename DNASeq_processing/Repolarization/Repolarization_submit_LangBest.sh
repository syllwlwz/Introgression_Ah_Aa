#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Repol_LangBest
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#./software/R-3.5.3/bin/Rscript Repolarization_Lang.R
grep "#" /prj/pflaphy-cutolgs/Filtered_lyrata/LangnewFil2.vcf | cat - /prj/pflaphy-cutolgs/GS_LangBest_fil/Langnewfil2_repol.vcf > /prj/pflaphy-cutolgs/GS_LangBest_fil/Langnewfil2_repolarized.vcf

#./software/R-3.5.3/bin/Rscript Repolarization_Best.R
#grep "#" /prj/pflaphy-cutolgs/Filtered_lyrata/BestnewFil2.vcf | cat - /prj/pflaphy-cutolgs/GS_LangBest_fil/Bestnewfil2_repol.vcf > /prj/pflaphy-cutolgs/GS_LangBest_fil/Bestnewfil2_repolarized.vcf
