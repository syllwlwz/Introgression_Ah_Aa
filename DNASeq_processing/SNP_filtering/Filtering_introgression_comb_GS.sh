#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_GS
#$ -l vf=65G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#MiasZapa
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 2 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasZapaarenosa.err

#Separate again
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/MiasZapaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Zapa_arenosa.args -O Filtered_introgression2/ZapaarenosaMFil.vcf\
# 2> Filtered_introgression2/ZapaarenosaMFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/MiasZapaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Mias_arenosa.args -O Filtered_introgression2/MiasarenosaZFil.vcf\
# 2> Filtered_introgression2/MiasarenosaZFil.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/ZapaarenosaMFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/ZapaarenosaM.table 2> Filtered_introgression2/ZapaarenosaMtable.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/MiasarenosaZFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/MiasarenosaZ.table 2> Filtered_introgression2/MiasarenosaZtable.err

#Separate populations
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_SNPs.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_SNPs.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_SNPs.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_SNPs.err

#Miashalleri
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 6 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/Miasarha1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Miasarha1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/Miasarha1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/Miasarha2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/Miasarha1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/Miasarha2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/Miasarha3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Miasarha3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Miasarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Miasarha.err

#Separate again
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Miasarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_halleri.args -O Filtered_introgression2/haMFil.vcf\
# 2> Filtered_introgression2/haMFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Miasarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Mias_arenosa.args -O Filtered_introgression2/MiasarenosahaFil.vcf\
# 2> Filtered_introgression2/MiasarenosahaFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haMFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/haMFil2.vcf\
# 2> Filtered_introgression2/haMFil2.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/MiasarenosahaFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/MiasarenosahaFil2.vcf\
# 2> Filtered_introgression2/MiasarenosahaFil2.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haMFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/MiasarenosahaFil2.vcf\
# -O Filtered_introgression2/haMFil3.vcf 2> Filtered_introgression2/haMFil3.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/MiasarenosahaFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/haMFil2.vcf -O Filtered_introgression2/MiasarenosahaFil3.vcf\
# 2> Filtered_introgression2/MiasarenosahaFil3.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/haMFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/haM.table 2> Filtered_introgression2/haMtable.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/MiasarenosahaFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/Miasarenosaha.table 2> Filtered_introgression2/Miasarenosahatable.err


#MiasKowa
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -R Alyrata_107.fa -env -minN 2 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/MiasKowaarenosa.err

#Separate again
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/MiasKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/KowaarenosaMFil.vcf\
# 2> Filtered_introgression2/KowaarenosaMFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/MiasKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Mias_arenosa.args -O Filtered_introgression2/MiasarenosaKFil.vcf\
# 2> Filtered_introgression2/MiasarenosaKFil.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KowaarenosaMFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/KowaarenosaM.table 2> Filtered_introgression2/KowaarenosaMtable.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/MiasarenosaKFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/MiasarenosaK.table 2> Filtered_introgression2/MiasarenosaKtable.err

#PiekKowa
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -R Alyrata_107.fa -env -minN 2 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/PiekKowaarenosa.err

#Separate again
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/PiekKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/KowaarenosaPFil.vcf\
# 2> Filtered_introgression2/KowaarenosaPFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/PiekKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Piek_arenosa.args -O Filtered_introgression2/PiekarenosaKFil.vcf\
# 2> Filtered_introgression2/PiekarenosaKFil.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KowaarenosaPFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/KowaarenosaP.table 2> Filtered_introgression2/KowaarenosaPtable.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/PiekarenosaKFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/PiekarenosaK.table 2> Filtered_introgression2/PiekarenosaKtable.err

#BukoKowa
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -R Alyrata_107.fa -env -minN 2 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/BukoKowaarenosa.err

#Separate again
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/BukoKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/KowaarenosaBFil.vcf\
# 2> Filtered_introgression2/KowaarenosaBFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/BukoKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Buko_arenosa.args -O Filtered_introgression2/BukoarenosaKFil.vcf\
# 2> Filtered_introgression2/BukoarenosaKFil.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KowaarenosaBFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/KowaarenosaB.table 2> Filtered_introgression2/KowaarenosaBtable.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/BukoarenosaKFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/BukoarenosaK.table 2> Filtered_introgression2/BukoarenosaKtable.err

#KatoKowa
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -R Alyrata_107.fa -env -minN 2 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/KatoKowaarenosa.err

#Separate again
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/KatoKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/KowaarenosaKFil.vcf\
# 2> Filtered_introgression2/KowaarenosaKFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/KatoKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kato_arenosa.args -O Filtered_introgression2/KatoarenosaKFil.vcf\
# 2> Filtered_introgression2/KatoarenosaKFil.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KowaarenosaKFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/KowaarenosaK.table 2> Filtered_introgression2/KowaarenosaKtable.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KatoarenosaKFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/KatoarenosaK.table 2> Filtered_introgression2/KatoarenosaKtable.err

#ZapaKowa
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -R Alyrata_107.fa -env -minN 2 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/ZapaKowaarenosa.err

#Separate again
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/ZapaKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/KowaarenosaZFil.vcf\
# 2> Filtered_introgression2/KowaarenosaZFil.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/ZapaKowaarenosa.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Zapa_arenosa.args -O Filtered_introgression2/ZapaarenosaKFil.vcf\
# 2> Filtered_introgression2/ZapaarenosaKFil.err

#Extract tables
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KowaarenosaZFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/KowaarenosaZ.table 2> Filtered_introgression2/KowaarenosaZtable.err

#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/ZapaarenosaKFil.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
# -O Filtered_introgression2/ZapaarenosaK.table 2> Filtered_introgression2/ZapaarenosaKtable.err



#Piekhalleri
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 6 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/Piekarha1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Piekarha1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/Piekarha1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/Piekarha2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/Piekarha1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/Piekarha2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/Piekarha3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Piekarha3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Piekarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Piekarha.err

#Separate again
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Piekarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_halleri.args -O Filtered_introgression2/haPFil.vcf\
 2> Filtered_introgression2/haPFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Piekarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Piek_arenosa.args -O Filtered_introgression2/PiekarenosahaFil.vcf\
 2> Filtered_introgression2/PiekarenosahaFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haPFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/haPFil2.vcf\
 2> Filtered_introgression2/haPFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/PiekarenosahaFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/PiekarenosahaFil2.vcf\
 2> Filtered_introgression2/PiekarenosahaFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haPFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/PiekarenosahaFil2.vcf\
 -O Filtered_introgression2/haPFil3.vcf 2> Filtered_introgression2/haPFil3.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/PiekarenosahaFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/haPFil2.vcf -O Filtered_introgression2/PiekarenosahaFil3.vcf\
 2> Filtered_introgression2/PiekarenosahaFil3.err

#Extract tables
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/haPFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/haP.table 2> Filtered_introgression2/haPtable.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/PiekarenosahaFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/Piekarenosaha.table 2> Filtered_introgression2/Piekarenosahatable.err



#Katohalleri
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kato_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 6 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/Katoarha1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Katoarha1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/Katoarha1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/Katoarha2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/Katoarha1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/Katoarha2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/Katoarha3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Katoarha3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Katoarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Katoarha.err

#Separate again
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Katoarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_halleri.args -O Filtered_introgression2/haKaFil.vcf\
 2> Filtered_introgression2/haKaFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Katoarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kato_arenosa.args -O Filtered_introgression2/KatoarenosahaFil.vcf\
 2> Filtered_introgression2/KatoarenosahaFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haKaFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/haKaFil2.vcf\
 2> Filtered_introgression2/haKaFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/KatoarenosahaFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/KatoarenosahaFil2.vcf\
 2> Filtered_introgression2/KatoarenosahaFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haKaFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/KatoarenosahaFil2.vcf\
 -O Filtered_introgression2/haKaFil3.vcf 2> Filtered_introgression2/haKaFil3.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/KatoarenosahaFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/haKaFil2.vcf -O Filtered_introgression2/KatoarenosahaFil3.vcf\
 2> Filtered_introgression2/KatoarenosahaFil3.err

#Extract tables
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/haKaFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/haKa.table 2> Filtered_introgression2/haKatable.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KatoarenosahaFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/Katoarenosaha.table 2> Filtered_introgression2/Katoarenosahatable.err



#Bukohalleri
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 6 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Bukoarha.err

#Separate again
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Bukoarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_halleri.args -O Filtered_introgression2/haBFil.vcf\
 2> Filtered_introgression2/haBFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Bukoarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Buko_arenosa.args -O Filtered_introgression2/BukoarenosahaFil.vcf\
 2> Filtered_introgression2/BukoarenosahaFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haBFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/haBFil2.vcf\
 2> Filtered_introgression2/haBFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/BukoarenosahaFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/BukoarenosahaFil2.vcf\
 2> Filtered_introgression2/BukoarenosahaFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haBFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/BukoarenosahaFil2.vcf\
 -O Filtered_introgression2/haBFil3.vcf 2> Filtered_introgression2/haBFil3.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/BukoarenosahaFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/haBFil2.vcf -O Filtered_introgression2/BukoarenosahaFil3.vcf\
 2> Filtered_introgression2/BukoarenosahaFil3.err

#Extract tables
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/haBFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/haB.table 2> Filtered_introgression2/haBtable.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/BukoarenosahaFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/Bukoarenosaha.table 2> Filtered_introgression2/Bukoarenosahatable.err



#Kowahalleri
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_arenosa_renamed_sorted_unphased_fil5_higherDP_20perc.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 6 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Kowaarha.err

#Separate again
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Kowaarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_halleri.args -O Filtered_introgression2/haKoFil.vcf\
 2> Filtered_introgression2/haKoFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Kowaarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Kowa_arenosa.args -O Filtered_introgression2/KowaarenosahaFil.vcf\
 2> Filtered_introgression2/KowaarenosahaFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haKoFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/haKoFil2.vcf\
 2> Filtered_introgression2/haKoFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/KowaarenosahaFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/KowaarenosahaFil2.vcf\
 2> Filtered_introgression2/KowaarenosahaFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haKoFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/KowaarenosahaFil2.vcf\
 -O Filtered_introgression2/haKoFil3.vcf 2> Filtered_introgression2/haKoFil3.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/KowaarenosahaFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/haKoFil2.vcf -O Filtered_introgression2/KowaarenosahaFil3.vcf\
 2> Filtered_introgression2/KowaarenosahaFil3.err

#Extract tables
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/haKoFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/haKo.table 2> Filtered_introgression2/haKotable.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/KowaarenosahaFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/Kowaarenosaha.table 2> Filtered_introgression2/Kowaarenosahatable.err

#Zapahalleri
#/prj/pflaphy-cutolgs/software/jre1.8.0_291/bin/java -Xmx50g -jar /prj/pflaphy-cutolgs/software/GenomeAnalysisTK37.jar -T CombineVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_arenosa_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Buko_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Mias_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Piek_halleri_renamed_sorted_unphased_fil5_higherDP.vcf -V /prj/pflaphy-gscan/Filtered_introgression2/Kowa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -V /prj/pflaphy-gscan/Filtered_introgression2/Zapa_halleri_renamed_sorted_unphased_fil5_higherDP.vcf\
# -R Alyrata_107.fa -env -minN 6 -filteredAreUncalled -o /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha1.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha1.err

#grep -v "#" /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha1.vcf | awk '{seen[$1,"\t",$2]++;lines[$1,"\t",$2]=$0}END{for(i in seen){if(seen[i]==1){print lines[i]}}}' | sort -k 1,1 -k2,2n > /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha2.vcf
#grep "#" /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha1.vcf | cat - /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha2.vcf > /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha3.vcf
#java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha3.vcf -R Alyrata_107.fa --select-type-to-include SNP -restrict-alleles-to BIALLELIC\
# --exclude-filtered -select "AF<1.0 && AF>0.0" --exclude-non-variants -O /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha.vcf 2> /prj/pflaphy-gscan/Filtered_introgression2/Zapaarha.err

#Separate again
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Zapaarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_halleri.args -O Filtered_introgression2/haZFil.vcf\
 2> Filtered_introgression2/haZFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/Zapaarha.vcf -R Alyrata_107.fa -sn Introgression_filter_list_Zapa_arenosa.args -O Filtered_introgression2/ZapaarenosahaFil.vcf\
 2> Filtered_introgression2/ZapaarenosahaFil.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haZFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/haZFil2.vcf\
 2> Filtered_introgression2/haZFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/ZapaarenosahaFil.vcf -R Alyrata_107.fa -select "AN>0" -O Filtered_introgression2/ZapaarenosahaFil2.vcf\
 2> Filtered_introgression2/ZapaarenosahaFil2.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/haZFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/ZapaarenosahaFil2.vcf\
 -O Filtered_introgression2/haZFil3.vcf 2> Filtered_introgression2/haZFil3.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -V Filtered_introgression2/ZapaarenosahaFil2.vcf -R Alyrata_107.fa -conc Filtered_introgression2/haZFil2.vcf -O Filtered_introgression2/ZapaarenosahaFil3.vcf\
 2> Filtered_introgression2/ZapaarenosahaFil3.err

#Extract tables
java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/haZFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/haZ.table 2> Filtered_introgression2/haZtable.err

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V Filtered_introgression2/ZapaarenosahaFil3.vcf -R Alyrata_107.fa -F CHROM -F POS -F AC -F AN -raw\
 -O Filtered_introgression2/Zapaarenosaha.table 2> Filtered_introgression2/Zapaarenosahatable.err



