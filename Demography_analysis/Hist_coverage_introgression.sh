#!/bin/bash

#bedtools version 2.26
#cohorts separately

#run from pflaphy-gscan

#$ -t 1-6
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Coverage
#$ -l vf=260G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat GROM.list | head -n $SGE_TASK_ID | tail -n 1)
./software/bedtools2/bin/bedtools coverage -hist -g Alyrata_107.fa.fai -a Alyrata_genes_sorted.gff -b GROM_analysis/$file'_arenosa.sorted.bam' -nonamecheck > Coverage/Hist_arenosa_$file.bed 2> Coverage/Hist_arenosa_$file.err
./software/bedtools2/bin/bedtools coverage -hist -g Alyrata_107.fa.fai -a Alyrata_genes_sorted.gff -b GROM_analysis/$file'_halleri.sorted.bam' -nonamecheck > Coverage/Hist_halleri_$file.bed 2> Coverage/Hist_halleri_$file.err
