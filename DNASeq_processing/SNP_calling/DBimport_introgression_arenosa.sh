#!/bin/bash

#run from pflaphy-gscan

#$ -t 1-9
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N DBimport
#$ -l vf=25G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

Scaffold=$(head -n $SGE_TASK_ID Mainscaffolds.list | tail -n 1)

java -Xmx80g -Xms80g -XX:ConcGCThreads=4 -jar software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport -R /prj/pflaphy-gscan/Alyrata_107.fa\
 -V /prj/pflaphy-gscan/HC1_introgression/Buko_a21.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Buko_a22.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Buko_a28.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Buko_a29.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Buko_a30.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Buko_a31.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Buko_a34b.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Buko_a34.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a21.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a22.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a23.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a24.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a27.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a28.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a29.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a30.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a33.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_a35.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h21.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h22.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h23.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h24.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h26.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h27.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h33.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kato_h35.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a04.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a05.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a06.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a07.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a08.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a09.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a11.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Kowa001a12.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias001a18.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias003a09.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias003a10.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias003a11.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias003a13.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias003a15.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias003a16.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_19a03x21a01_d.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_23a03x22a03_k.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a12xa05_l.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a42.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a43.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a44.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a45.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a46.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a47.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a48.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a50.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Mias_a55.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a10.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a12.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a13.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a14.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a15.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a16.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a17.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a18.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a19.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a1.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a20-1.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a20-2.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a22.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a2.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a4.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a5.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a7.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a8.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_a9.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_h14_1.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Piek_h6_2.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa002a01.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa002a02.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa002a03.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa002a04.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa002a05.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa002a07.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa004a11.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa008a19.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa_11a03x12a01_l.g.vcf.gz -V /prj/pflaphy-gscan/HC1_introgression/Zapa_a09xa03_b.g.vcf.gz\
 --use-jdk-inflater --tmp-dir /prj/pflaphy-gscan/HC2_temp/ -DF NotDuplicateReadFilter -L $Scaffold --genomicsdb-workspace-path /prj/pflaphy-gscan/DB_introgression_arenosa_$Scaffold/ --batch-size 50 --max-num-intervals-to-import-in-parallel 4 --merge-input-intervals true --genomicsdb-shared-posixfs-optimizations true 2> /prj/pflaphy-gscan/HC2_introgression/GenomicsDBimport_arenosa_$Scaffold.err


