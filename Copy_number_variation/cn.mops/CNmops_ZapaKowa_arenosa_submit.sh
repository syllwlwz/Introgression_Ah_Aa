#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cnmops
#$ -l vf=12G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

/prj/pflaphy-gscan/software/R-4.1.0/bin/Rscript CNmops_introgression_Zapa_arenosa.R 2> Cnmops_ZapaKowaarenosa.err
