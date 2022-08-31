#!/bin/bash

#$ -t 1-60
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N RNASeqqc_lyrata_errorcor
#$ -l vf=7G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Cutadapt_introrna.list | head -n $SGE_TASK_ID | tail -n 1)
./software/qualimap_v2.2.1/qualimap rnaseq -a proportional -bam Error_corrected2/$file'_outfile.bam' -gtf Alyrata_384_v2.1.gene_exons.gtf -p strand-specific-reverse -pe 2> Error_corrected2/$file.errorcor.rnaseqqc.err

