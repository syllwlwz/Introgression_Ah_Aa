#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Merge
#$ -l vf=56G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#./software/samtools-1.11/samtools merge -f /prj/pflaphy-gscan/dedup/Bko_h1_7.dedup.bam /prj/pflaphy-gscan/dedup/Buko_h01.dedup.bam /prj/pflaphy-gscan/dedup/Buko_h07.dedup.bam 2> /prj/pflaphy-gscan/dedup/Buko_h1_7comb.err
#./software/samtools-1.11/samtools merge -f /prj/pflaphy-gscan/dedup/Zako_h01comb.dedup.bam /prj/pflaphy-gscan/dedup/Zako_h01.dedup.bam /prj/pflaphy-gscan/dedup/Zako002h01.dedup.bam 2> /prj/pflaphy-gscan/dedup/Zako_h01comb.err
#./software/samtools-1.11/samtools merge -f /prj/pflaphy-gscan/dedup/Zako_h09comb.dedup.bam /prj/pflaphy-gscan/dedup/Zako_h09.dedup.bam /prj/pflaphy-gscan/dedup/Zako001h01.dedup.bam 2> /prj/pflaphy-gscan/dedup/Zako_h09comb.err
#./software/samtools-1.11/samtools merge -f /prj/pflaphy-gscan/dedup/Kowa_h04comb.dedup.bam /prj/pflaphy-gscan/dedup/Kowa_h04.dedup.bam /prj/pflaphy-gscan/dedup/Kowa003h04.dedup.bam 2> /prj/pflaphy-gscan/dedup/Kowa_h04comb.err
#./software/samtools-1.11/samtools merge -f /prj/pflaphy-gscan/dedup/Kowa_h02comb.dedup.bam /prj/pflaphy-gscan/dedup/Kowa_h02.dedup.bam /prj/pflaphy-gscan/dedup/Kowa003h02.dedup.bam 2> /prj/pflaphy-gscan/dedup/Kowa_h02comb.err
#./software/samtools-1.11/samtools merge -f /prj/pflaphy-gscan/dedup/Zapa_h11comb.dedup.bam /prj/pflaphy-gscan/dedup/Zapa_h11.dedup.bam /prj/pflaphy-gscan/dedup/Zapa001h11.dedup.bam 2> /prj/pflaphy-gscan/dedup/Zapa_h11comb.err
#./software/samtools-1.11/samtools merge -f /prj/pflaphy-gscan/dedup/Zapa_h01comb.dedup.bam /prj/pflaphy-gscan/dedup/Zapa_h01.dedup.bam /prj/pflaphy-gscan/dedup/Zapa001h01.dedup.bam 2> /prj/pflaphy-gscan/dedup/Zapa_h01comb.err

#java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Zako_h01comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Zako_h01comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Zako_h01 RGPL=illumina RGPU=AAAAAA RGSM=Zako_h01\
# 2> /prj/pflaphy-gscan/dedup/Zako_h01comb.namefixed.dedup.err

#java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Zako_h09comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Zako_h09comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Zako_h09 RGPL=illumina RGPU=AAAAAA RGSM=Zako_h09\ 
# 2> /prj/pflaphy-gscan/dedup/Zako_h09comb.namefixed.dedup.err

#java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Kowa_h04comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Kowa_h04comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Kowa_h04 RGPL=illumina RGPU=AAAAAA RGSM=Kowa_h04\ 
# 2> /prj/pflaphy-gscan/dedup/Kowa_h04comb.namefixed.dedup.err

#java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Kowa_h02comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Kowa_h02comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Kowa_h02 RGPL=illumina RGPU=AAAAAA RGSM=Kowa_h02\ 
# 2> /prj/pflaphy-gscan/dedup/Kowa_h02comb.namefixed.dedup.err

#java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Zapa_h11comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Zapa_h11comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Zapa_h11 RGPL=illumina RGPU=AAAAAA RGSM=Zapa_h11\ 
# 2> /prj/pflaphy-gscan/dedup/Zapa_h11comb.namefixed.dedup.err

#java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Zapa_h01comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Zapa_h01comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Zapa_h01 RGPL=illumina RGPU=AAAAAA RGSM=Zapa_h01\ 
# 2> /prj/pflaphy-gscan/dedup/Zapa_h01comb.namefixed.dedup.err


#java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Buko_h1_7.dedup.bam O=/prj/pflaphy-gscan/dedup/Buko_h1_7.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Buko_h01 RGPL=illumina RGPU=AAAAAA RGSM=Buko_h01\ 
# 2> /prj/pflaphy-gscan/dedup/Buko_h1_7.namefixed.dedup.err

java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Piek_h05comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Piek_h05comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Piek_h05 RGPL=illumina RGPU=AAAAAA RGSM=Piek_h05\ 
 2> /prj/pflaphy-gscan/dedup/Piek_h05comb.namefixed.dedup.err


java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=/prj/pflaphy-gscan/dedup/Piek_h06comb.dedup.bam O=/prj/pflaphy-gscan/dedup/Piek_h06comb.namefixed.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Piek_h06comb RGPL=illumina RGPU=AAAAAA RGSM=Piek_h06comb\ 
 2> /prj/pflaphy-gscan/dedup/Piek_h06comb.namefixed.dedup.err
