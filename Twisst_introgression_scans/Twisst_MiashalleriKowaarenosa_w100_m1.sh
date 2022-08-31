#!/bin/bash
#Twisst
#runs per scaffold

#$ -t 1-8
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Twisst
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

i=$(cat Mainscaffolds.list | head -n $SGE_TASK_ID | tail -n 1)
python3 software/twisst/twisst.py -t Phyml_lyrata/MiashalleriKowaarenosalyr.Vcf.scf$i.w100m1.LDphase.phyml_bionj.trees.gz -w Twisst/MiashalleriKowaarenosa.weights.w100.m1.$i.csv.gz --outputTopos Twisst/MiashalleriKowaarenosa.topologies.w100.m1.$i.trees\
 --outgroup A_ly\
 -g A_ly -g Mias_arenosa -g Halleri -g Kowa_arenosa --method complete --groupsFile Twisst/MiashalleriKowaarenosa_Twisst_groups2.txt --distsFile Twisst/MiashalleriKowaarenosa.w100.m1.$i.dist 2> Twisst/MiashalleriKowaarenosa_Twisst_w100_m1_$i.err

