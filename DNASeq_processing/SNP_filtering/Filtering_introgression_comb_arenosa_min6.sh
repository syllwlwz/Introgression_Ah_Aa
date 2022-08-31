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
#DPmax=2*$mode
#DPmin=0.5*$mode

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil2.vcf -R Alyrata_107.fa\
# --filter-name "DP" -filter "DP<"$DPmin" || DP>"$DPmax"" -O /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.err

#separate into populations for AN filter
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Buko_arenosa.args -O Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil4.vcf 2> Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kato_arenosa.args -O Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil4.vcf 2> Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Mias_arenosa.args -O Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil4.vcf 2> Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Piek_arenosa.args -O Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil4.vcf 2> Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Zapa_arenosa.args -O Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil4.vcf 2> Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil4.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_arenosa_renamed_sorted_unphased_fil3.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil4.vcf 2> Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil4.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5.vcf\
# 2> Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5.vcf\
# 2> Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5.vcf\
# 2> Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5.vcf\
# 2> Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.vcf\
# 2> Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=32" -V Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5.vcf\
# 2> Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5.err

Buko_n=$(wc -l Introgression_filter_list_Buko_arenosa.args | cut -f 1)
Mias_n=$(wc -l Introgression_filter_list_Mias_arenosa.args | cut -f 1)
Piek_n=$(wc -l Introgression_filter_list_Piek_arenosa.args | cut -f 1)
Zapa_n=$(wc -l Introgression_filter_list_Zapa_arenosa.args | cut -f 1)
Kowa_n=$(wc -l Introgression_filter_list_Kowa_arenosa.args | cut -f 1)
Kato_n=$(wc -l Introgression_filter_list_Kato_arenosa.args | cut -f 1)

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=24" -V Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_min6.vcf\
 2> Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_min6.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=24" -V Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_min6.vcf\
 2> Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_min6.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=24" -V Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_min6.vcf\
 2> Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_min6.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=24" -V Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_min6.vcf\
 2> Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_min6.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -select "AN>=24" -V Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil4.vcf -O Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_min6.vcf\
 2> Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_min6.err









#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr.err

#cat Introgression_filter_list_Buko_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Zapa_arenosa.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_BukoZapa.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_BukoZapa.args -O /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr.vcf -R Alyrata_107.fa -conc /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa.vcf\
# -O /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr_fil.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr_fil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.05" --exclude-non-variants -sn Introgression_filter_list_BukoZapa.args -O /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa2.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa2.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr.vcf -R Alyrata_107.fa -conc /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa2.vcf\
# -O /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr_fil2.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr_fil2.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr1.err

#cat Introgression_filter_list_Mias_arenosa.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Zapa_arenosa.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_MiasZapa.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiasZapa.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiasZapa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasZapa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiasZapa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr_fil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr_fil.vcf -R Alyrata_107.fa\
# --exclude-filtered -select "AF<1.0 && AF>0.05" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr_fil_MAF.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasZapalyr_fil_MAF.err


#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/PiekZapalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/PiekZapalyr1.err

#cat Introgression_filter_list_Piek_arenosa.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Zapa_arenosa.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_PiekZapa.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/PiekZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_PiekZapa.args -O /prj/pflaphy-gscan/Filtered_introgression2/PiekZapa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/PiekZapa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/PiekZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/PiekZapa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/PiekZapalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/PiekZapalyr_fil.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/KatoZapalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/KatoZapalyr1.err

#cat Introgression_filter_list_Kato_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Zapa_arenosa.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_KatoZapa.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/KatoZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_KatoZapa.args -O /prj/pflaphy-gscan/Filtered_introgression2/KatoZapa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/KatoZapa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/KatoZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/KatoZapa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/KatoZapalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/KatoZapalyr_fil.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.err

#cat Introgression_filter_list_Buko_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Zapa_arenosa.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_BukoZapa.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_BukoZapa.args -O /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa.err#

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/BukoZapa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoZapalyr_fil.err






#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasKowalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasKowalyr1.err

#cat Introgression_filter_list_Mias_arenosa.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Kowa_arenosa.args Introgression_filter_list_Kowa_halleri.args > Introgression_filter_list_MiasKowa.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasKowalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiasKowa.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiasKowa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasKowa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasKowalyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiasKowa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiasKowalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasKowalyr_fil.err

