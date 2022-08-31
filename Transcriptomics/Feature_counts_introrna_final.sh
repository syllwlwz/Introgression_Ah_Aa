#!/bin/bash

#$ -t 1-60
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Feature_counts
#$ -l vf=200M
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

file=$(cat Cutadapt_introrna.list | head -n $SGE_TASK_ID | tail -n 1)
./software/subread-2.0.3-Linux-x86_64/bin/featureCounts -a Alyrata_384_v2.1.gene_exons.gtf -o Counts_FC/$file.counts.txt -F 'GTF' -t exon -g gene_id -O -M --fraction -p --countReadPairs -T 12 -s 2 --fracOverlap 0.15 --largestOverlap Error_corrected2/$file'_outfile.bam' 2> Counts_FC/$file.counts.err
#-F fiel format of annotation
#-t count over exons
#-g summarize counts over gene_ids
#-O all overlapping meta-features
#each read is counted once even if exons are present several times due to splicing variants
#--minOverlap minimum overlapping bp of read for counting counted for both reads
#--fracOverlap minimum fraction of read overlapping feature counted for both reads
#-M count multiple mapping reads
#--fraction "proportional" counting for multiple mapping reads
#-p read pairs
#--countReadPairs count fragments, if only 1 read of pair mapped also counted
#-T number of threads
#-s 2 strand-specific reverse
