#!/bin/bash

#cutadapt version 2.10

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N LDprune
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

python3 from ./software/ScanTools import scantools; project = scantools("/prj/software/Scantools/Popkey_introgression.csv", "encoding"); project.splitVCFs(vcf_dir=Structure, min_dp=0, mffg=1); 
./software/LD_prune.py -v Structure -w 500 -r 10 -d 1000 -o Structure -gz false -s true -c 0.1 -m 1 -mf 0 
#parser.add_argument('-v', type = str, metavar = 'vcf_path', required = True, help = 'path to vcfs')
#parser.add_argument('-w', type = str, metavar = 'Window_Size', required = True, help = 'size of scaffold window')
#parser.add_argument('-r', type = str, metavar = 'Number_Replications', required = True, help = 'Number of replicate data sets')
#parser.add_argument('-d', type = str, metavar = 'Window_Distance', required = True, help = 'distance between any 2 windows')
#parser.add_argument('-m', type = str, metavar = 'Missing_Data', required = True, help = 'amount of missing data to allow per site (between 0-1)')
#parser.add_argument('-mf', type = str, metavar = 'Minimum_Frequency', required = True, help = 'minimum minor allele frequency')
#parser.add_argument('-o', type = str, metavar = 'Output_Prefix', required = True, help = 'Vcfs retain original scaffold name but the concatenated Structure input file will be a text file with specified by output and within the LD_pruned directory')
#parser.add_argument('-s', type = str, metavar = 'Subset?', required = False, default = 'false', help = 'if true, this will subsample polyploid data to create psuedo-diploid data')
#parser.add_argument('-c', type = float, metavar = 'maximum_correlation', required = True, default = '1.0', help = 'Maximum correlation between adjacent sites allowed for site to be used')
#parser.add_argument('-gz', type = str, metavar = 'gzipped?', required = True, help = 'are vcfs gzipped (true) or not (false)')


