#!/bin/bash

#run from pflaphy-gscan

#$ -t 1-204
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N nQuire_depth_histo
#$ -l vf=500M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Total_samples_introgression.list | head -n $SGE_TASK_ID | tail -n 1)

./software/nQuire/nQuire histo nQuire_depth/$file'_denoised'.bin > nQuire_depth/$file.histo.table

