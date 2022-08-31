#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Repolarization_key
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#python3 software/ScanTools/repolarization_lookupKey_introgression.py -v /prj/pflaphy-gscan/Repolarized/Outgroup_unphased_bi.vcf -o /prj/pflaphy-gscan/Repolarized/Introgression -mi 4 -mp 0.7 2> /prj/pflaphy-gscan/Repolarized/Repolarization_key.err
#python3 software/ScanTools/repolarization_lookupKey_introgression.py -v /prj/pflaphy-gscan/Repolarized/Outgroup_unphased_bi.vcf -o /prj/pflaphy-gscan/Repolarized/Introgression_05 -mi 4 -mp 0.5 2> /prj/pflaphy-gscan/Repolarized/Repolarization_key_05.err
python3 software/ScanTools/repolarization_lookupKey_introgression.py -v /prj/pflaphy-gscan/Repolarized/Outgroup_unphased_bi.vcf -o /prj/pflaphy-gscan/Repolarized/Introgression_03 -mi 4 -mp 0.3 2> /prj/pflaphy-gscan/Repolarized/Repolarization_key_03.err

python3 software/ScanTools/repolarization_lookupKey_introgression.py -v /prj/pflaphy-gscan/Repolarized/Outgroup_unphased_bi.vcf -o /prj/pflaphy-gscan/Repolarized/Introgression_mi5 -mi 5 -mp 0.7 2> /prj/pflaphy-gscan/Repolarized/Repolarization_key_mi5.err
python3 software/ScanTools/repolarization_lookupKey_introgression.py -v /prj/pflaphy-gscan/Repolarized/Outgroup_unphased_bi.vcf -o /prj/pflaphy-gscan/Repolarized/Introgression_mi5_05 -mi 5 -mp 0.5 2> /prj/pflaphy-gscan/Repolarized/Repolarization_key_05_mi5.err
python3 software/ScanTools/repolarization_lookupKey_introgression.py -v /prj/pflaphy-gscan/Repolarized/Outgroup_unphased_bi.vcf -o /prj/pflaphy-gscan/Repolarized/Introgression_mi5_03 -mi 5 -mp 0.3 2> /prj/pflaphy-gscan/Repolarized/Repolarization_key_03_mi5.err
