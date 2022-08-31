#!/bin/bash

#Error Correction
#script_comex

#$ -t 56-60
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N errorCorrection
#$ -l vf=45G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Cutadapt_introrna.list | head -n $SGE_TASK_ID | tail -n 1)

InputSam=mapped_lyrata2/$file.sam
if [ -f $InputSam ];
	then
		echo "File $InputSam exists."
	else
		echo "File $InputSam does not exist."
		exit
fi

IFS=\.

arr=($InputSam)


input=$file
echo $input

IFS=$OIFS

referenceGenome=Alyrata_107.fa

if [ -f $referenceGenome ];
	then
		echo "File $referenceGenome exists."
	else
		echo "File $referenceGenome does not exist."
		exit
fi

python2 ./toprintend1.py mapped_lyrata2/$input'.sam' ./Error_corrected2/$input'_end.sam'

echo "End created"

python2 ./Selectnonrepeated1.py ./Error_corrected2/$input'_end.sam' ./Error_corrected2/$input'_hits_nonrepeated.sam'

echo "non-repeated hits finished"

python2 ./removeEnd1.py ./Error_corrected2/$input'_hits_nonrepeated.sam' >./Error_corrected2/$input'_output_final.sam'

echo "End removed"

./software/samtools-1.14/samtools view -bT $referenceGenome ./Error_corrected2/$input'_output_final.sam' >./Error_corrected2/$input'_outfile.bam'

echo "All tophat-errors corrected for $input."

rm ./Error_corrected2/$input'_end.sam'
rm ./Error_corrected2/$input'_hits_nonrepeated.sam'
rm ./Error_corrected2/$input'_output_final.sam'
echo "Temporary files have been deleted"


