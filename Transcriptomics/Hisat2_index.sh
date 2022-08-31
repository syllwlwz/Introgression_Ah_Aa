#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Hisat2_index
#$ -l vf=10G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 6

./software/hisat2-2.2.0/hisat2-build -f Alyrata_107.fa Alyrata_all.index --ss Splice_sites_Alyrata_all.hisat2 --exon Exons_Alyrata_all.hisat2 -p 6 2> Hisat2_index.err
