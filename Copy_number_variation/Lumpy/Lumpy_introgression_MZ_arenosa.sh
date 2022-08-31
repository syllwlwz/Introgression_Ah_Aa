#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Lumpy
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

Input1=Mias_arenosa
Input2=Zapa_arenosa



RL=$(software/samtools-1.11/samtools view dedup/$file.dedup.bam | awk '{print length($10)}' | head -1000 | sort -nu | tail -n 1)

./software/samtools-1.11/samtools view dedup/$file.dedup.bam | tail -n+100000 | software/lumpy-sv/scripts/pairend_distro.py -r $RL -X 4 -N 10000 -o dedup/$file.lib1.histo > dedup/$file.distrib
lumpy \
    -mw 4 \
    -tt 0 \
    -pe id:sample1,bam_file:sample1.discordants.bam,read_group:rg1,read_group:rg2,histo_file:sample1.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -pe id:sample1,bam_file:sample1.discordants.bam,read_group:rg3,histo_file:sample1.lib2.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -pe id:sample2,bam_file:sample2.discordants.bam,read_group:rg4,histo_file:sample2.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:sample1,bam_file:sample1.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:sample2,bam_file:sample2.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
    > multi_sample.vcf
