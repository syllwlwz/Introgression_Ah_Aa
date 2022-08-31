#!/bin/bash

#Deduplication
#picard tools 2.21.9

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Cat_structure_vcfs
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -Xmx50g -jar /prj/pflaphy-gscan/software/picard.jar GatherVcfs O=/prj/pflaphy-gscan/Structure/LD_Pruned/All5_introgression_4dgsites.LD_Pruned.vcf I=All5_4dg_scaffold_1.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_2.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_3.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_4.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_5.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_6.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_7.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_8.rep0.LD_Pruned.vcf I=All5_4dg_scaffold_9.rep0.LD_Pruned.vcf 2> /prj/pflaphy-gscan/Structure/LD_Pruned/Gathervcfs2.err
