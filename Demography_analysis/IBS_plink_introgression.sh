#!/bin/bash


#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N plink_IBS
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

./software/plink --distance bin 1-ibs --bfile Structure/LD_Pruned/All_introgression_4dgsites.LD_Pruned.renamed.plink --out Structure/LD_Pruned/All_introgression_4dgsites.LD_Pruned.renamed.plink.IBS --allow-extra-chr 2> Structure/LD_Pruned/IBS_plink.err

