#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_allintrogression
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#Filter 1: select bialleleic SNPs without missing genotypes
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V HC2_introgression/Introgression_renamed_sorted.vcf.gz -R Alyrata_107.fa --select-type-to-include SNP --exclude-non-variants -restrict-alleles-to BIALLELIC -O Filtered_introgression/All_introgressionFil1_pops.vcf 2> Filtered_introgression/All_introgressionFil1_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V Filtered_introgression/All_introgressionFil1_pops.vcf -R Alyrata_107.fa -G-filter "DP<4.0" -G-filter-name "lowDP" -O Filtered_introgression/All_introgressionFil2_pops.vcf 2> Filtered_introgression/All_introgressionFil2_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil2_pops.vcf -R Alyrata_107.fa --set-filtered-gt-to-nocall -O Filtered_introgression/All_introgressionFil3_pops.vcf 2> Filtered_introgression/All_introgressionFil3_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V Filtered_introgression/All_introgressionFil3_pops.vcf -R Alyrata_107.fa --filter-name "DPgeno" -filter "QD<2.0" --filter-name "FS" -filter "FS>40.0" --filter-name "MQ" -filter "MQ<50.0" --filter-name "MQRanksum" -filter "MQRankSum < -2.5" --filter-name "ReadosRankSum" -filter "ReadPosRankSum < - 4.0" --filter-name "SOR" -filter "SOR>4.0" -O Filtered_introgression/All_introgressionFil5a_pops.vcf 2> Filtered_introgression/All_introgressionFil5a_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil5a_pops.vcf -R Alyrata_107.fa --exclude-filtered -O Filtered_introgression/All_introgressionFil5_pops.vcf 2> Filtered_introgression/All_introgressionFil5_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil5_pops.vcf -R Alyrata_107.fa -F DP -O Filtered_introgression/All_introgressionFil5.table 2> Filtered_introgression/All_introgressionFil5_DPtable_pops.err

#grep -v "DP" Filtered_introgression/All_introgressionFil5.table > Filtered_introgression/All_introgressionFil5_DP_hQD.table
#mode=$(./software/R-4.0.4/bin/Rscript -e 'library("LaplacesDemon"); data<-read.table("Filtered_introgression/All_introgressionFil5_DP_hQD.table",header=F); cat(Mode(data$V1)[1])')
#DPmax=2*$mode
#DPmin=0.5*$mode

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -V Filtered_introgression/All_introgressionFil5_pops.vcf -R Alyrata_107.fa --filter-name "DP" -filter "DP<"$DPmin" || DP>"$DPmax"" -O Filtered_introgression/All_introgressionFil6_pops.vcf 2> Filtered_introgression/All_introgressionFil6_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil6_pops.vcf -R Alyrata_107.fa --exclude-filtered -O Filtered_introgression/All_introgressionFil7_pops.vcf 2> Filtered_introgression/All_introgressionFil7_pops.err

#separate into pops
#select % missingness per pop
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Buko.args -O Filtered_introgression/All_introgressionFil8_Buko.vcf 2> Filtered_introgression/All_introgressionFil8_Buko.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kato.args -O Filtered_introgression/All_introgressionFil8_Kato.vcf 2> Filtered_introgression/All_introgressionFil8_Kato.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa.args -O Filtered_introgression/All_introgressionFil8_Kowa.vcf 2> Filtered_introgression/All_introgressionFil8_Kowa.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Mias.args -O Filtered_introgression/All_introgressionFil8_Mias.vcf 2> Filtered_introgression/All_introgressionFil8_Mias.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Piek.args -O Filtered_introgression/All_introgressionFil8_Piek.vcf 2> Filtered_introgression/All_introgressionFil8_Piek.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_outgroup.args -O Filtered_introgression/All_introgressionFil8_outgroup.vcf 2> Filtered_introgression/All_introgressionFil8_outgroup.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Zapa.args -O Filtered_introgression/All_introgressionFil8_Zapa.vcf 2> Filtered_introgression/All_introgressionFil8_Zapa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_Buko.vcf -R Alyrata_107.fa --exclude-filtered -select "AN > 20" -O Filtered_introgression/All_introgressionFil9_Buko.vcf 2> Filtered_introgression/All_introgressionFil9_Buko.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_Kato.vcf -R Alyrata_107.fa --exclude-filtered -select "AN > 20" -O Filtered_introgression/All_introgressionFil9_Kato.vcf 2> Filtered_introgression/All_introgressionFil9_Kato.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_Kowa.vcf -R Alyrata_107.fa --exclude-filtered -select "AN > 20" -O Filtered_introgression/All_introgressionFil9_Kowa.vcf 2> Filtered_introgression/All_introgressionFil9_Kowa.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_Mias.vcf -R Alyrata_107.fa --exclude-filtered -select "AN > 20" -O Filtered_introgression/All_introgressionFil9_Mias.vcf 2> Filtered_introgression/All_introgressionFil9_Mias.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_Piek.vcf -R Alyrata_107.fa --exclude-filtered -select "AN > 20" -O Filtered_introgression/All_introgressionFil9_Piek.vcf 2> Filtered_introgression/All_introgressionFil9_Piek.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_outgroup.vcf -R Alyrata_107.fa --exclude-filtered -select "AN > 20" -O Filtered_introgression/All_introgressionFil9_outgroup.vcf 2> Filtered_introgression/All_introgressionFil9_outgroup.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_Zapa.vcf -R Alyrata_107.fa --exclude-filtered -select "AN > 20" -O Filtered_introgression/All_introgressionFil9_Zapa.vcf 2> Filtered_introgression/All_introgressionFil9_Zapa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar IndexFeatureFile -I Filtered_introgression/All_introgressionFil9_Buko.vcf

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil7_pops.vcf -R Alyrata_107.fa -conc Filtered_introgression/All_introgressionFil9_Buko.vcf --exclude-filtered -O Filtered_introgression/All_introgressionFil8_pops.vcf 2> Filtered_introgression/All_introgressionFil8_pops.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil8_pops.vcf -R Alyrata_107.fa -conc Filtered_introgression/All_introgressionFil9_Kato.vcf --exclude-filtered -O Filtered_introgression/All_introgressionFil9_pops.vcf 2> Filtered_introgression/All_introgressionFil9_pops.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil9_pops.vcf -R Alyrata_107.fa -conc Filtered_introgression/All_introgressionFil9_Kowa.vcf --exclude-filtered -O Filtered_introgression/All_introgressionFil10_pops.vcf 2> Filtered_introgression/All_introgressionFil10_pops.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil10_pops.vcf -R Alyrata_107.fa -conc Filtered_introgression/All_introgressionFil9_Mias.vcf --exclude-filtered -O Filtered_introgression/All_introgressionFil11_pops.vcf 2> Filtered_introgression/All_introgressionFil11_pops.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil11_pops.vcf -R Alyrata_107.fa -conc Filtered_introgression/All_introgressionFil9_Piek.vcf --exclude-filtered -O Filtered_introgression/All_introgressionFil12_pops.vcf 2> Filtered_introgression/All_introgressionFil12_pops.err
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil12_pops.vcf -R Alyrata_107.fa -conc Filtered_introgression/All_introgressionFil9_Zapa.vcf --exclude-filtered -O Filtered_introgression/All_introgressionFil13_pops.vcf 2> Filtered_introgression/All_introgressionFil13_pops.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil13_pops.vcf -R Alyrata_107.fa --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O Filtered_introgression/All_introgressionFil14_pops.vcf 2> Filtered_introgression/All_introgressionFil14_pops.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Filtered_introgression/All_introgressionFil_pops.table 2> Filtered_introgression/All_introgressionFiltable_pops.err

#split into pops for structure
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Buko.args -O Filtered_introgression/All_introgressionFil14_pops_Buko.vcf 2> Filtered_introgression/All_introgressionFil14_pops_Buko.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kato.args -O Filtered_introgression/All_introgressionFil14_pops_Kato.vcf 2> Filtered_introgression/All_introgressionFil14_pops_Kato.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa.args -O Filtered_introgression/All_introgressionFil14_pops_Kowa.vcf 2> Filtered_introgression/All_introgressionFil14_pops_Kowa.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Mias.args -O Filtered_introgression/All_introgressionFil14_pops_Mias.vcf 2> Filtered_introgression/All_introgressionFil14_pops_Mias.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Piek.args -O Filtered_introgression/All_introgressionFil14_pops_Piek.vcf 2> Filtered_introgression/All_introgressionFil14_pops_Piek.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Zapa.args -O Filtered_introgression/All_introgressionFil14_pops_Zapa.vcf 2> Filtered_introgression/All_introgressionFil14_pops_Zapa.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression/All_introgressionFil14_pops.vcf -R Alyrata_107.fa -sn Introgression_filter_list_outgroup.args -O Filtered_introgression/All_introgressionFil14_pops_outgroup.vcf 2> Filtered_introgression/All_introgressionFil14_pops_outgroup.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops_Buko.vcf -R Alyrata_107.fa -F CHROM -F POS -F REF -F AN -F DP -GF GT -O Filtered_introgression/All_introgressionFil14_pops_Buko_raw.table 2> Filtered_introgression/All_introgressionFil14_pops_Buko_raw_table.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops_Kato.vcf -R Alyrata_107.fa -F CHROM -F POS -F REF -F AN -F DP -GF GT -O Filtered_introgression/All_introgressionFil14_pops_Kato_raw.table 2> Filtered_introgression/All_introgressionFil14_pops_Kato_raw_table.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops_Kowa.vcf -R Alyrata_107.fa -F CHROM -F POS -F REF -F AN -F DP -GF GT -O Filtered_introgression/All_introgressionFil14_pops_Kowa_raw.table 2> Filtered_introgression/All_introgressionFil14_pops_Kowa_raw_table.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops_Mias.vcf -R Alyrata_107.fa -F CHROM -F POS -F REF -F AN -F DP -GF GT -O Filtered_introgression/All_introgressionFil14_pops_Mias_raw.table 2> Filtered_introgression/All_introgressionFil14_pops_Mias_raw_table.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops_Piek.vcf -R Alyrata_107.fa -F CHROM -F POS -F REF -F AN -F DP -GF GT -O Filtered_introgression/All_introgressionFil14_pops_Piek_raw.table 2> Filtered_introgression/All_introgressionFil14_pops_Piek_raw_table.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops_Zapa.vcf -R Alyrata_107.fa -F CHROM -F POS -F REF -F AN -F DP -GF GT -O Filtered_introgression/All_introgressionFil14_pops_Zapa_raw.table 2> Filtered_introgression/All_introgressionFil14_pops_Zapa_raw_table.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression/All_introgressionFil14_pops_outgroup.vcf -R Alyrata_107.fa -F CHROM -F POS -F REF -F AN -F DP -GF GT -O Filtered_introgression/All_introgressionFil14_pops_outgroup_raw.table 2> Filtered_introgression/All_introgressionFil14_pops_outgroup_raw_table.err




