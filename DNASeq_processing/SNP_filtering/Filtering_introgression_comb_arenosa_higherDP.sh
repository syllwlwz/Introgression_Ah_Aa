#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_arenosa
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#Annotate genotypes with depth below 6 as missing and variants failing the best practices
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased.vcf -R Alyrata_107.fa -G-filter "DP<6.0" -G-filter-name "lowDP"\
# --set-filtered-genotype-to-no-call --filter-name "DPgeno" -filter "QD<2.0" --filter-name "FS" -filter "FS>40.0" --filter-name "MQ" -filter "MQ<50.0" --filter-name "MQRanksum" -filter "MQRankSum < -2.5" --filter-name "ReadosRankSum" -filter "ReadPosRankSum < - 4.0"\
# --filter-name "SOR" -filter "SOR>4.0" -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil1.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil1.vcf -R Alyrata_107.fa\
# --exclude-filtered -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.err

#Annotate loci with too low/ excessive read depth less/ more than 2x the mode of the depth distribution
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.vcf -R Alyrata_107.fa -F DP\
# -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.table 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.table.err

#grep -v "DP" /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.table | grep -v "NA" | grep -v "0" > /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2_DP.table
#mode=$(./software/R-4.1.0/bin/Rscript -e 'library("LaplacesDemon"); data<-read.table("/prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2_DP.table",header=F); cat(Mode(data$V1[data$V1>100])[1])')
#DPmax=8*$mode
#DPmin=0.5*$mode

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.vcf -R Alyrata_107.fa\
# --filter-name "DP" -filter "DP>"$DPmax"" -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.err

#separate into populations for AN filter
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Buko_arenosa.args -O Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil4_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kato_arenosa.args -O Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil4_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Mias_arenosa.args -O Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil4_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Piek_arenosa.args -O Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil4_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Zapa_arenosa.args -O Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil4_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil4_higherDP.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# 2> Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# 2> Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# 2> Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# 2> Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# 2> Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# 2> Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Mias_SNPs.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Mias_SNPs.err

