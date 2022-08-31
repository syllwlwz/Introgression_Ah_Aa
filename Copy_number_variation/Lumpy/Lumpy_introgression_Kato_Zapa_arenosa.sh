#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Lumpy
#$ -l vf=30G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

means=()
stdevs=()
files=()
pe=()
sr=()
rd=()
for file in $(cat dedup/Kato_arenosa.list)
do
files+=($file)
while IFS=':' read -ra line
do
    l=($(echo ${line[@]} | tr ":" "\t"))
    mean=${l[1]}
#    means+=($mean)
    stdev=${l[3]}
#    stdevs+=($stdev)
done < dedup/$file.distrib

#echo ${means[@]}
#echo ${stdevs[@]}
#echo ${files[@]}
rd=$(grep $file[[:blank:]] dedup/Read_lengths_samples.txt | cut -f 2)

pe+=(-pe id:$file,bam_file:dedup/$file.discordants,histo_file:dedup/$file.lib1.histo,mean:$mean,stdev:$stdev,read_length:$rd,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:25 \ )
sr+=(-sr id:$file,bam_file:dedup/$file.splitters,back_distance:10,weight:1,min_mapping_threshold:25 \ )

done

echo ${pe[@]}
echo ${sr[@]}

means_Zapa=()
stdevs_Zapa=()
files_Zapa=()
pe_Zapa=()
sr_Zapa=()
rd_Zapa=()
for file_Zapa in $(cat dedup/Zapa_arenosa.list)
do
files_Zapa+=($file_Zapa)
while IFS=':' read -ra line
do
    l_Zapa=($(echo ${line[@]} | tr ":" "\t"))
    mean_Zapa=${l[1]}
#    means+=($mean)
    stdev_Zapa=${l[3]}
#    stdevs+=($stdev)
done < dedup/$file_Zapa.distrib

#echo ${means[@]}
#echo ${stdevs[@]}
#echo ${files[@]}
rd_Zapa=$(grep $file_Zapa[[:blank:]] dedup/Read_lengths_samples.txt | cut -f 2)

pe_Zapa+=(-pe id:$file_Zapa,bam_file:dedup/$file_Zapa.discordants,histo_file:dedup/$file_Zapa.lib1.histo,mean:$mean_Zapa,stdev:$stdev_Zapa,read_length:$rd_Zapa,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:25 \ )
sr_Zapa+=(-sr id:$file_Zapa,bam_file:dedup/$file_Zapa.splitters,back_distance:10,weight:1,min_mapping_threshold:25 \ )

done

 ./software/lumpy-sv/bin/lumpy \
    -mw 4 \
    -tt 0 \
    ${pe[@]}\
    ${pe_Zapa[@]}\
    ${sr[@]}\
    ${sr_Zapa[@]}\
    > Lumpy/Kato_Zapa_arenosa.vcf 2> Lumpy/Kato_Zapa_arenosa.txt

