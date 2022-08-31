#!/bin/bash

#run from pflaphy-gscan

#$ -t 26
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N nQuire_depth
#$ -l vf=50M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 9

file=$(cat /prj/pflaphy-gscan/Outgroup_renamed.list | head -n $SGE_TASK_ID | tail -n 1)

max_dp=$(grep $file'.lyrata.mean_DOC_all.sample_statistics' /prj/pflaphy-gscan/dedup/Coverage_threshold.table | cut -f 2)
echo $max_dp
/prj/pflaphy-gscan/software/nQuire/nQuire create -b /prj/pflaphy-gscan/dedup/$file.dedup.bam -o /prj/pflaphy-gscan/nQuire_Ute2/$file -x -q 30 -c 15 -m $max_dp -r /prj/pflaphy-gscan/Mainscaffolds.bed -y > /prj/pflaphy-gscan/nQuire_Ute2/$file.out 2> /prj/pflaphy-gscan/nQuire_Ute2/$file.err
/prj/pflaphy-gscan/software/nQuire/nQuire denoise /prj/pflaphy-gscan/nQuire_Ute2/$file-bedcc.bin -o /prj/pflaphy-gscan/nQuire_Ute2/$file'_denoised' > /prj/pflaphy-gscan/nQuire_Ute2/$file.denoised.out 2> /prj/pflaphy-gscan/nQuire_Ute2/$file.denoised.err
/prj/pflaphy-gscan/software/nQuire/nQuire lrdmodel -t 9 /prj/pflaphy-gscan/nQuire_Ute2/$file'_denoised'.bin > /prj/pflaphy-gscan/nQuire_Ute2/$file.lrdmodel.table
#higher delta log-likelihood of a fixed model = lower support for the corresponding ploidy level

