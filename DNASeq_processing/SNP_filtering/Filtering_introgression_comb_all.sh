#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Combine_all
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#########################################################
#Whole dataset
#########################################################
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_thaliana_renamed_sorted_unphased_fil4.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_croatica_renamed_sorted_unphased_fil4.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_cebennensis_renamed_sorted_unphased_fil4.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_pedemontana_renamed_sorted_unphased_fil4.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 16 -o /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/All2_woKatoha.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/All2_woKatoha.vcf > /prj/pflaphy-gscan/Filtered_introgression2/All3_woKatoha.vcf

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All3_woKatoha.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_halleri_renamed_sorted_unphased_fil4_higherDP.vcf -R Alyrata_107.fa\
# --exclude-filtered -conc /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/Kato_halleri2b.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Kato_hallerib.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants\
# -V /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_halleri2b.vcf -R Alyrata_107.fa -env -o /prj/pflaphy-gscan/Filtered_introgression2/All2.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All2.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/All3.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All3.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/All2.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/All2b.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/All2.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/All2b.vcf > /prj/pflaphy-gscan/Filtered_introgression2/All4.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All4.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/All5.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All5.err




#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants\
# -V /prj/pflaphy-gscan/Structure/VCF_storage/All5_4dg.vcf -V /prj/pflaphy-gscan/Structure/600_dataset_may2020_arenosa_snp_raw.snps.fourfold.dp8nc.m0.5.vcf.gz\
# -R Alyrata_107.fa -env -minN 2 -o /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb.vcf 2> /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb.err

#grep -v "#" /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_combb.vcf
#grep "#" /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb.vcf | cat - /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_combb.vcf > /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb4.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb4.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb5.vcf 2> /prj/pflaphy-gscan/Structure/All5_600arenosa_4dg_comb5.err



##########################################################
#Without Arabidopsis outgroup#
##########################################################

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
## -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 11 -o /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha_wooutgroup.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha_wooutgroup.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha_wooutgroup.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/All2_woKatoha_wooutgroup.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/All1_woKatoha_wooutgroup.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/All2_woKatoha_wooutgroup.vcf > /prj/pflaphy-gscan/Filtered_introgression2/All3_woKatoha_wooutgroup.vcf

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All3_woKatoha_wooutgroup.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha_wooutgroup.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha_wooutgroup.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_halleri_renamed_sorted_unphased_fil4_higherDP.vcf -R Alyrata_107.fa\
 --exclude-filtered -conc /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha_wooutgroup.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/Kato_halleri2b_wooutgroup.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Kato_hallerib_wooutgroup.err

/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants\
 -V /prj/pflaphy-gscan/Filtered_introgression2/All_woKatoha_wooutgroup.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_halleri2b_wooutgroup.vcf -R Alyrata_107.fa -env -o /prj/pflaphy-gscan/Filtered_introgression2/All2_wooutgroup.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All2_wooutgroup.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All2_wooutgroup.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
 --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/All3_wooutgroup.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All3_wooutgroup.err

grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/All2_wooutgroup.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/All2b_wooutgroup.vcf
grep "#" /prj/pflaphy-gscan/Filtered_introgression2/All2_wooutgroup.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/All2b_wooutgroup.vcf > /prj/pflaphy-gscan/Filtered_introgression2/All4_wooutgroup.vcf
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/All4_wooutgroup.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
 -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/All5_wooutgroup.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/All5_wooutgroup.err


