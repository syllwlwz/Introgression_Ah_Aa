#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_Twisst
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 4 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr2.>
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr3.vcf

#cat Introgression_filter_list_Mias_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Zapa_arenosa.args > Introgression_filter_list_MiasBukoZapaar.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiasBukoZapaar.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaar.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaar.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaar.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasBukoZapaarlyr_fil.err


#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 9 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/Miashal>
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr3.vcf

#cat Introgression_filter_list_Mias_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Kowa_arenosa.args\
# Introgression_filter_list_Kowa_halleri.args Introgression_filter_list_Zapa_halleri.args Introgression_filter_list_Zapa_arenosa.args > Introgression_filter_list_MiashalleriKowaZapaarenosa.args
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiashalleriKowaZapaarenosa.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Miashalle>

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaZapaarenosalyr_fil.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 8 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr3.vcf

#cat Introgression_filter_list_Mias_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Kowa_arenosa.args\
# Introgression_filter_list_Kowa_halleri.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_MiashalleriKowaarenosa.args

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiashalleriKowaarenosa.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosalyr_fil.err


#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 8 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k\
# 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr3.vcf

#cat Introgression_filter_list_Piek_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Kowa_arenosa.args\
# Introgression_filter_list_Kowa_halleri.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_PiekhalleriKowaarenosa.args

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_PiekhalleriKowaarenosa.args -O /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosa.vcf 2>\
# /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/PiekhalleriKowaarenosalyr_fil.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 8 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k\
# 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr3.vcf

#cat Introgression_filter_list_Buko_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Kowa_arenosa.args\
# Introgression_filter_list_Kowa_halleri.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_BukohalleriKowaarenosa.args

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_BukohalleriKowaarenosa.args -O /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosa.vcf 2>\
# /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/BukohalleriKowaarenosalyr_fil.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 8 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k\
# 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr3.vcf

#cat Introgression_filter_list_Kato_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Kowa_arenosa.args\
# Introgression_filter_list_Kowa_halleri.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_KatohalleriKowaarenosa.args

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_KatohalleriKowaarenosa.args -O /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosa.vcf 2>\
# /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosa.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/KatohalleriKowaarenosalyr_fil.err

/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
 -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
 -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
 -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
 -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
 -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
 -R Alyrata_107.fa -env -minN 8 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr1.err

grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k\
 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr2.vcf
grep "#" /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr3.vcf

cat Introgression_filter_list_Zapa_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Kowa_arenosa.args\
 Introgression_filter_list_Kowa_halleri.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_ZapahalleriKowaarenosa.args

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
 --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_ZapahalleriKowaarenosa.args -O /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosa.vcf 2>\
 /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosa.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
 --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosa.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr_fil.vcf\
 2> /prj/pflaphy-gscan/Filtered_introgression2/ZapahalleriKowaarenosalyr_fil.err





#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_thaliana_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 8 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath3.vcf

#cat Introgression_filter_list_Mias_arenosa.args Introgression_filter_list_Buko_halleri.args Introgression_filter_list_Mias_halleri.args Introgression_filter_list_Piek_halleri.args Introgression_filter_list_Kowa_arenosa.args\
# Introgression_filter_list_Kowa_halleri.args Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_MiashalleriKowaarenosa.args

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiashalleriKowaarenosa.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosa2.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosa2.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosa2.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriKowaarenosath_fil.err


#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriZapaarenosath_fil.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn A_th_2 -sn Mias_19a03x21a01_d -sn Mias009h03 -sn Zapa_11a03x12a01_l -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarth_set1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarth_set1.err


#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiashalleriZapaarenosath_fil.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Mias_19a03x21a01_d -sn Mias009h03 -sn Zapa_11a03x12a01_l -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaar_set1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaar_set1.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarth_set1.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaar_set1.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarth_set1_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarth_set1_fil.err



#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_lyrata_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr3.vcf

#cat Introgression_filter_list_Mias_arenosa.argsIntrogression_filter_list_Mias_halleri.args Introgression_filter_list_Kowa_arenosa.args\
# Introgression_filter_list_Kowa_halleri.args > Introgression_filter_list_MiasarhaKowaarha.args

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiasarhaKowaarha.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarha.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarha.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaKowaarhalyr_fil.err

#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Introgression_thaliana_renamed_sorted_unphased_fil4.vcf\
# -R Alyrata_107.fa -env -minN 5 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/Miasarhaowaarhath3.vcf

#cat Introgression_filter_list_Mias_arenosa.argsIntrogression_filter_list_Mias_halleri.args Introgression_filter_list_Zapa_arenosa.args\
# Introgression_filter_list_Zapa_halleri.args > Introgression_filter_list_MiasarhaZapaarha.args

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -sn Introgression_filter_list_MiasarhaZapaarha.args -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarha.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -conc /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarha.vcf -O /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath_fil.vcf\
# 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasarhaZapaarhath_fil.err




