#!/bin/bash

#$ -t 9
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Phyml
#$ -l vf=1G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

i=$(cat Mainscaffolds.list | head -n $SGE_TASK_ID | tail -n 1)
python3 software/genomics_general/phylo/phyml_sliding_windows2.py -g Twisst/MiashalleriKowaarenosalyrata_SNPs_fil2.geno.gz --include $i --windType sites -w 100 -M 1 --optimise n --maxLDphase --indFile Twisst/MiashalleriKowaarenosalyr.txt\
 -p Phyml_lyrata/MiashalleriKowaarenosalyr.Vcf.scf$i.w100m1.LDphase.phyml_bionj --phyml /prj/pflaphy-gscan/software/phyml-3.3.20190321/src/phyml --log Phyml_lyrata/Phyml_MiashalleriKowaarenosalyr_log_w100_$i.txt -T 12\
 --tmp Phyml_lyrata/Phyml_MiashalleriKowaarenosalyr_temp_w100_$i.txt 2> Phyml_lyrata/Phyml_MiashalleriKowaarenosalyr_w100_$i.err
