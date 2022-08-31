#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_allintrogression
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_croatica_renam
#ed_sorted_unphased_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_cebennensis_renamed_sorted_unphased_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_pedemontana_renamed_sorted_unphased_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_thali
#ana_renamed_sorted_unphased_higherDP.vcf\
# -R Alyrata_107.fa -env -o /prj/pflaphy-gscan/Filtered_introgression2/All_higherDP.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All.err

#Annotate genotypes with depth below 6 as missing and variants failing the best practices
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_higherDP.vcf -R Alyrata_107.fa -G-filter "DP<6.0" -G-filter-name "lowDP"\
# --set-filtered-genotype-to-no-call --filter-name "DPgeno" -filter "QD<2.0" --filter-name "FS" -filter "FS>40.0" --filter-name "MQ" -filter "MQ<50.0" --filter-name "MQRanksum" -filter "MQRankSum < -2.5" --filter-name "ReadosRankSum" -filter "ReadPosRankSum < - 4.0"\
# --filter-name "SOR" -filter "SOR>4.0" -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil1_higherDP.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil1.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil1_higherDP.vcf -R Alyrata_107.fa\
# --exclude-filtered -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2_higherDP.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2.err

#Annotate loci with too low/ excessive read depth less/ more than 2x the mode of the depth distribution
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2_higherDP.vcf -R Alyrata_107.fa -F DP\
# -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2.table 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2.table.err

#grep -v "DP" /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2.table | grep -v "NA" | grep -v "0" > /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2_DP.table
#mode=$(./software/R-4.1.0/bin/Rscript -e 'library("LaplacesDemon"); data<-read.table("/prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2_DP.table",header=F); cat(Mode(data$V1[data$V1>100])[1])')
#DPmax=8*$mode
#DPmin=0.5*$mode

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil2.vcf -R Alyrata_107.fa\
# --filter-name "DP" -filter "DP>"$DPmax"" -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil3_higherDP.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil3_higherDP.err

#separate into populations for AN filter
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Buko_halleri.args -O Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Mias_halleri.args -O Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Piek_halleri.args -O Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Zapa_halleri.args -O Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil4.err

Buko_n=$(wc -l Introgression_filter_list_Buko_halleri.args | cut -f 1 -d " ")
Mias_n=$(wc -l Introgression_filter_list_Mias_halleri.args | cut -f 1 -d " ")
Piek_n=$(wc -l Introgression_filter_list_Piek_halleri.args | cut -f 1 -d " ")
Zapa_n=$(wc -l Introgression_filter_list_Zapa_halleri.args | cut -f 1 -d " ")
Kowa_n=$(wc -l Introgression_filter_list_Kowa_halleri.args | cut -f 1 -d " ")

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=(0.8*$Buko_n*2)" -V Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
 2> Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=(0.8*$Mias_n*2)" -V Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
 2> Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=(0.8*$Piek_n*2)" -V Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
 2> Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=(0.8*$Zapa_n*2)" -V Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
 2> Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP_20perc.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_halleri_renamed_sorted_unphased_fil3_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kowa_halleri.args -O Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil4_higherDP.vcf 2> Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=16" -V Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil4_higherDP.vcf -O Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# 2> Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5.err








#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar IndexFeatureFile -I Filtered_introgression/Halleri_introgressionFil9_Zapa_higherDP.vcf


#merge species
#select biall. SNPs
#subset to sample subsets for Twisst with removing non-variants -env







#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V HC2_introgression/Introgression_halleri_renamed_sorted_higherDP.vcf.gz -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC -O Filtered_introgression/Halleri_introgressionFil1_pops_higherDP.vcf 2> Filtered_introgression/Halleri_introgressionFil1_pops.err

#select % missingness per pop
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/Halleri_introgressionFil7_pops_higherDP.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Buko_halleri.args -O Filtered_introgression/Halleri_introgressionFil8_Buko_higherDP.vcf 2> Filtered_introgression/Halleri_introgressionFil8_Buko.err

#./software/bcftools-1.12/bcftools view -i '(N_SAMPLES-N_MISSING)>7' Filtered_introgression/Halleri_introgressionFil8_Buko_higherDP.vcf > Filtered_introgression/Halleri_introgressionFil9_Buko_higherDP.vcf 2> Filtered_introgression/Halleri_introgressionFil9_Buko.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar IndexFeatureFile -I Filtered_introgression/Halleri_introgressionFil9_Zapa_higherDP.vcf

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/Halleri_introgressionFil7_pops_higherDP.vcf -R Alyrata_107.fa -conc Filtered_introgression/Halleri_introgressionFil9_Buko_higherDP.vcf --exclude-filtered -O Filtered_introgression/Halleri_introgressionFil8_pops_higherDP.vcf 2> Filtered_introgression/Halleri_introgressionFil8_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/Halleri_introgressionFil13_pops_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O Test/Halleri_introgressionFil14_pops_higherDP.vcf 2> Test/Halleri_introgressionFil14_pops.err
#./software/jdk-11.0.11/bin/java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/Halleri_introgressionFil13_pops_higherDP.vcf -R Alyrata_107.fa --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O Test/Halleri_introgressionFil14_pops2_higherDP.vcf 2> Test/Halleri_introgressionFil14_pops2.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/Halleri_introgressionFil14_pops_higherDP.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Filtered_introgression/Halleri_introgressionFil_pops.table 2> Filtered_introgression/Halleri_introgressionFiltable_pops.err
