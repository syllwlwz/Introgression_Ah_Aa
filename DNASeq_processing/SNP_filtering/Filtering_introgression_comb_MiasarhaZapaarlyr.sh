#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_MiasarhaZapaarlyr
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5.vcf\
 -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5.vcf\
 -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
 -R Alyrata_107.fa -env -minN 4 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarlyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarlyr1.err

cat Introgression_filter_list_Mias_arenosa.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Zapa_arenosa.args > Introgression_filter_list_MiasarhaZapaar.args
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarlyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
 --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiasarhaZapaar.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaar.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaar.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarlyr1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
 --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaar.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarlyr_fil.vcf\
 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarlyr_fil.err

