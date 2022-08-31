#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N STRUCTURE
#$ -l vf=200M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 8

/homes/lsyllwas/.local/bin/structure_threader run -R 3 -o /prj/pflaphy-gscan/Structure/ -st /homes/lsyllwas/.local/bin/structure -K 15 -i /prj/pflaphy-gscan/Structure/LD_Pruned/LD_struc.StructureInput.rep1.LD_Pruned.Diploidized.txt --params /homes/lsyllwas/.local/bin/mainparams -t 8 2> /prj/pflaphy-gscan/Structure/Structure.err
# --use-ind-labels
