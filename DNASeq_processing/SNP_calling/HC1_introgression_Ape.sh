##!/bin/bash

#$ -t 1-2
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N HC1_introgression_Ape
#$ -l vf=14G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

file=$(cat A_pe.args | head -n $SGE_TASK_ID | tail -n 1)
echo $file

java -Xmx50g -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar HaplotypeCaller -I /prj/pflaphy-gscan/dedup/$file.dedup.bam --min-base-quality-score 25 --minimum-mapping-quality 25 -DF NotDuplicateReadFilter -R /prj/pflaphy-gscan/Alyrata_107.fa -O /prj/pflaphy-gscan/HC1_introgression/$file.g.vcf.gz -ploidy 2 --pcr-indel-model NONE --heterozygosity 0.005 --indel-heterozygosity 0.001 --emit-ref-confidence GVCF --output-mode EMIT_ALL_ACTIVE_SITES --add-output-vcf-command-line 2> /prj/pflaphy-gscan/HC1_introgression/$file.err

# min-base quality score raufgesetzt auf 25(DFEFAULT 10)
# mimimum mapping quality ebenfalls auf 25(DEFAULT 20)
# -df NotDuplicateReadFilter -> ausschalten von Readfilter der Duplicates entfernt
#-RF GoodCigarReadFilter, notsecondaryalignmentReadFilter -> eigentlich obsolet da er auch automatisch angewandt wird (Vgl GATK seite )
#pcr indel model -> None da PCR free libraries
