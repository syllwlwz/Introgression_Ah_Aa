#!/bin/bash

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Lumpy
#$ -l vf=50G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

means=()
stdevs=()
files=()
pe=()
sr=()
rd=()
for file in $(cat dedup/Buko_halleri.list) 
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

means_Kowa=()
stdevs_Kowa=()
files_Kowa=()
pe_Kowa=()
sr_Kowa=()
rd_Kowa=()
for file_Kowa in $(cat dedup/Kowa_halleri.list) 
do
files_Kowa+=($file_Kowa)
while IFS=':' read -ra line
do
    l_Kowa=($(echo ${line[@]} | tr ":" "\t"))
    mean_Kowa=${l[1]}
#    means+=($mean)
    stdev_Kowa=${l[3]}
#    stdevs+=($stdev)
done < dedup/$file_Kowa.distrib

#echo ${means[@]}
#echo ${stdevs[@]}
#echo ${files[@]}
rd_Kowa=$(grep $file_Kowa[[:blank:]] dedup/Read_lengths_samples.txt | cut -f 2)

pe_Kowa+=(-pe id:$file_Kowa,bam_file:dedup/$file_Kowa.discordants,histo_file:dedup/$file_Kowa.lib1.histo,mean:$mean_Kowa,stdev:$stdev_Kowa,read_length:$rd_Kowa,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:25 \ )
sr_Kowa+=(-sr id:$file_Kowa,bam_file:dedup/$file_Kowa.splitters,back_distance:10,weight:1,min_mapping_threshold:25 \ )

done

 ./software/lumpy-sv/bin/lumpy \
    -mw 4 \
    -tt 0 \
    ${pe[@]}\
    ${pe_Kowa[@]}\
    ${sr[@]}\
    ${sr_Kowa[@]}\
    > Lumpy/Buko_Kowa_halleri.vcf 2> Lumpy/Buko_Kowa_halleri.txt

