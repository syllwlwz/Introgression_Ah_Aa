#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Pop_sep_repol
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#separate into populations and species
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Buko_halleri.args\
 -O Repolarized/Buko_halleri_all.vcf 2> Repolarized/Buko_halleri_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Mias_halleri.args\
 -O Repolarized/Mias_halleri_all.vcf 2> Repolarized/Mias_halleri_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Piek_halleri.args\
 -O Repolarized/Piek_halleri_all.vcf 2> Repolarized/Piek_halleri_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Zapa_halleri.args\
 -O Repolarized/Zapa_halleri_all.vcf 2> Repolarized/Zapa_halleri_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Kato_h09\
 -O Repolarized/Kato_halleri_all.vcf 2> Repolarized/Kato_halleri_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kowa_halleri.args\
 -O Repolarized/Kowa_halleri_all.vcf 2> Repolarized/Kowa_halleri_repol_all.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Buko_arenosa.args\
 -O Repolarized/Buko_arenosa_all.vcf 2> Repolarized/Buko_arenosa_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Mias_arenosa.args\
 -O Repolarized/Mias_arenosa_all.vcf 2> Repolarized/Mias_arenosa_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Piek_arenosa.args\
 -O Repolarized/Piek_arenosa_all.vcf 2> Repolarized/Piek_arenosa_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Zapa_arenosa.args\
 -O Repolarized/Zapa_arenosa_all.vcf 2> Repolarized/Zapa_arenosa_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kato_arenosa.args\
 -O Repolarized/Kato_arenosa_all.vcf 2> Repolarized/Kato_arenosa_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_Kowa_arenosa.args\
 -O Repolarized/Kowa_arenosa_all.vcf 2> Repolarized/Kowa_arenosa_repol_all.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_thaliana.args\
 -O Repolarized/thaliana_all.vcf 2> Repolarized/thaliana_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/thaliana_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/thaliana_repol_all.vcf 2> Repolarized/thaliana_repol2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_lyrata.args\
 -O Repolarized/lyrata_all.vcf 2> Repolarized/lyrata_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/lyrata_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/lyrata_repol_all.vcf 2> Repolarized/lyrata_repol2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_croatica.args\
 -O Repolarized/croatica_all.vcf 2> Repolarized/croatica_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/croatica_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/croatica_repol_all.vcf 2> Repolarized/croatica_repol2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_pedemontana.args\
 -O Repolarized/pedemontana_all.vcf 2> Repolarized/pedemontana_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/pedemontana_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/pedemontana_repol_all.vcf 2> Repolarized/pedemontana_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_repolarized.vcf -R Alyrata_107.fa --exclude-filtered -sn Introgression_filter_list_cebennensis.args\
 -O Repolarized/cebennensis_all.vcf 2> Repolarized/cebennensis_repol_all.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/cebennensis_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/cebennensis_all_Treemix.vcf 2> Repolarized/cebennensis_all_repol2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Buko_halleri_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Buko_halleri_repol_all.vcf 2> Repolarized/Buko_halleri_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Mias_halleri_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Mias_halleri_repol_all.vcf 2> Repolarized/Mias_halleri_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Piek_halleri_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Piek_halleri_repol_all.vcf 2> Repolarized/Piek_halleri_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Zapa_halleri_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Zapa_halleri_repol_all.vcf 2> Repolarized/Zapa_halleri_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Kato_halleri_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Kato_halleri_repol_all.vcf 2> Repolarized/Kato_halleri_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Kowa_halleri_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Kowa_halleri_repol_all.vcf 2> Repolarized/Kowa_halleri_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Buko_arenosa_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Buko_arenosa_repol_all.vcf 2> Repolarized/Buko_arenosa_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Mias_arenosa_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Mias_arenosa_repol_all.vcf 2> Repolarized/Mias_arenosa_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Piek_arenosa_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Piek_arenosa_repol_all.vcf 2> Repolarized/Piek_arenosa_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Zapa_arenosa_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Zapa_arenosa_repol_all.vcf 2> Repolarized/Zapa_arenosa_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Kato_arenosa_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Kato_arenosa_repol_all.vcf 2> Repolarized/Kato_arenosa_repol2.err
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Repolarized/Kowa_arenosa_all.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw -O Repolarized/Kowa_arenosa_repol_all.vcf 2> Repolarized/Kowa_arenosa_repol2.err

