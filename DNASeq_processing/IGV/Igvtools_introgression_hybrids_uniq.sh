##!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N igvtools_introgression
#$ -l vf=8G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

./software/IGV_2.11.2/igvtools count --minMapQuality 0 -w 100 --includeDuplicates aligned/Mias001a17.bowtie2.arhacomb.uniq.sort.bam aligned/Mias001a17.bowtie2.arhacomb.uniq.sort.bam.MQ0.tdf /prj/pflaphy-gscan/Halleri_arenosa.genome\
 2> ./software/IGV_2.11.2/Tdf.MQ0.Mias001a17.bowtie2.arhacomb.uniq.err
./software/IGV_2.11.2/igvtools count --minMapQuality 0 -w 100 --includeDuplicates aligned/Mias001a22.bowtie2.arhacomb.uniq.sort.bam aligned/Mias001a22.bowtie2.arhacomb.uniq.sort.bam.MQ0.tdf /prj/pflaphy-gscan/Halleri_arenosa.genome\
 2> ./software/IGV_2.11.2/Tdf.MQ0.Mias001a22.bowtie2.arhacomb.uniq.err
./software/IGV_2.11.2/igvtools count --minMapQuality 0 -w 100 --includeDuplicates aligned/Zapa002a06.bowtie2.arhacomb.uniq.sort.bam aligned/Zapa002a06.bowtie2.arhacomb.uniq.sort.bam.MQ0.tdf /prj/pflaphy-gscan/Halleri_arenosa.genome\
 2> ./software/IGV_2.11.2/Tdf.MQ0.Zapa002a06.bowtie2.arhacomb.uniq.err
