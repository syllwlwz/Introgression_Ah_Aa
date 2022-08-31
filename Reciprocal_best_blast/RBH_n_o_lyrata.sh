#!/bin/bash

#reciprocal best blast hits

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N RBH_n_o
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

./software/rbh-master/rbh -1 Alyrata_384_v2.1.cds_primaryTranscriptOnly.fa -2 TAIR10_cds_20110103_representative_gene_model_updated.txt --no-tblastx -t 12 -o --outfmt6 -b software/ncbi-blast-2.9.0+/bin -r RBH_n_o_lyrata 2> RBH_n_o_lyrata/RBH_lyrata.err

