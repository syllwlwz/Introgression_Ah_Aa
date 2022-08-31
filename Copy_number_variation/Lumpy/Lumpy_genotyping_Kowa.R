filelist = list.files(pattern = "*Kowa_arenosa.vcf$") # Make a file list from all count text files

for (file in filelist)
	{assign(file,read.table(file,sep="\t",header=F))
	}

for (file in filelist)
	{file_now<-get(file)
	colnames(file_now)[1:9]=c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names")
	file_now1<-file_now[grep("SVTYPE=DUP",file_now$INFO),]
	file_now2<-file_now[grep("SVTYPE=DEL",file_now$INFO),]
	file_now3<-file_now[grep("SVTYPE=CNV",file_now$INFO),]
	file_now4<-file_now[grep("SVTYPE=INS",file_now$INFO),]

	assign(paste0(file,"_DUP"),file_now1)
	assign(paste0(file,"_DEL"),file_now2)
	assign(paste0(file,"_CNV"),file_now3)
	assign(paste0(file,"_INS"),file_now4)
	}

##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">
##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description="Secondary breakend in a multi-line variants">
##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">
##INFO=<ID=BD,Number=.,Type=Integer,Description="Amount of BED evidence supporting the variant across all samples">
##INFO=<ID=EV,Number=.,Type=String,Description="Type of LUMPY evidence contributing to the variant call">
##INFO=<ID=PRPOS,Number=.,Type=String,Description="LUMPY probability curve of the POS breakend">
##INFO=<ID=PREND,Number=.,Type=String,Description="LUMPY probability curve of the END breakend">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
##FORMAT=<ID=BD,Number=1,Type=Integer,Description="Amount of BED evidence supporting the variant">


#from Lumpy forum: In Delly, SVs are flagged as "PASS" if >=3 paired-ends support the variant.
#You can filter variants for quality by requiring SU>=threshold.

#split INFO column to quality filter
require(tidyr)

Lists<-list.files(pattern = "*.list")
for (list in Lists)
	{assign(list,scan(list,quote="",what="character"))
	}

require(vegan)
for (file in filelist)
	{file_now<-get(paste0(file,"_DUP"))
	test<-separate(data=file_now,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="left")
	test$SU<-gsub("SU=","",test$NA10)
	test$PE<-gsub("PE=","",test$NA11)
	test$SR<-gsub("SR=","",test$NA12)

	list_now<-get(gsub(".vcf",".list",(gsub("_Kowa_","_",file))))
	species<-gsub(".vcf","",strsplit(file,"_")[[1]][3])
	Kowa_list<-get(paste0("Kowa_",species,".list"))
	colnames(file_now)[c(10:length(colnames(file_now)))]<-c(list_now,Kowa_list)

	SU_M<-rep(0,nrow(file_now))
	SR_M<-rep(0,nrow(file_now))
	PE_M<-rep(0,nrow(file_now))
	SU_NM<-rep(0,nrow(file_now))
	SR_NM<-rep(0,nrow(file_now))
	PE_NM<-rep(0,nrow(file_now))

#cumulative evidence in target population
	for (column in 10:(10+length(list_now)-1))
		{file_now2<-file_now
		column_now<-separate(data=file_now2,into=c("GT","SU","PE","SR"),col=column,sep=":")
		SU_M<-SU_M+as.numeric(as.character(column_now$SU))
		SR_M<-SR_M+as.numeric(as.character(column_now$SR))
		PE_M<-PE_M+as.numeric(as.character(column_now$PE))
		}
	for (column in (10+length(list_now)):(10-1+length(list_now)+length(Kowa_list)))
		{file_now2<-file_now
		column_now<-separate(data=file_now2,into=c("GT","SU","PE","SR"),col=column,sep=":")
		SU_NM<-SU_NM+as.numeric(as.character(column_now$SU))
		SR_NM<-SR_NM+as.numeric(as.character(column_now$SR))
		PE_NM<-PE_NM+as.numeric(as.character(column_now$PE))
		}

	file_now$PE<-test$PE
	file_now$SR<-test$SR
	file_now$SU<-test$SU
	file_now$PE_M<-PE_M
	file_now$SU_M<-SU_M
	file_now$SR_M<-SR_M
	file_now$PE_NM<-PE_NM
	file_now$SU_NM<-SU_NM
	file_now$SR_NM<-SR_NM

	assign(paste0(file,"_DUP"),file_now)


#	pdf("Lumpy_parameters_Qfil_logtenplusone.pdf",width=10,height=10,paper="special")
#	cloud(log10(test2$PE+1) ~ log10(test2$SR+1)*log10(test2$SU+1))
#	cloud(log10(test2$SU+1) ~ log10(test2$SR+1)*log10(test2$PE+1))
#	cloud(log10(test2$SR+1) ~ log10(test2$PE+1)*log10(test2$SU+1))
#	dev.off()

#summary(log10(test2$SR+1))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.5799  1.1139  3.9028 
#summary(log10(test2$PE+1))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   0.699   1.079   1.112   1.580   3.893 
#summary(log10(test2$SU+1))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6990  0.9542  1.3010  1.3891  1.7634  3.9570

	}

for (file in filelist)
	{file_now<-get(paste0(file,"_DEL"))
	test<-separate(data=file_now,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="left")
	test$SU<-gsub("SU=","",test$NA10)
	test$PE<-gsub("PE=","",test$NA11)
	test$SR<-gsub("SR=","",test$NA12)

	list_now<-get(gsub(".vcf",".list",(gsub("_Kowa_","_",file))))
	species<-gsub(".vcf","",strsplit(file,"_")[[1]][3])
	Kowa_list<-get(paste0("Kowa_",species,".list"))
	colnames(file_now)[c(10:length(colnames(file_now)))]<-c(list_now,Kowa_list)

	SU_M<-rep(0,nrow(file_now))
	SR_M<-rep(0,nrow(file_now))
	PE_M<-rep(0,nrow(file_now))
	SU_NM<-rep(0,nrow(file_now))
	SR_NM<-rep(0,nrow(file_now))
	PE_NM<-rep(0,nrow(file_now))

#cumulative evidence in target population
	for (column in 10:(10+length(list_now)-1))
		{file_now2<-file_now
		column_now<-separate(data=file_now2,into=c("GT","SU","PE","SR"),col=column,sep=":")
		SU_M<-SU_M+as.numeric(as.character(column_now$SU))
		SR_M<-SR_M+as.numeric(as.character(column_now$SR))
		PE_M<-PE_M+as.numeric(as.character(column_now$PE))
		}
	for (column in (10+length(list_now)):(10-1+length(list_now)+length(Kowa_list)))
		{file_now2<-file_now
		column_now<-separate(data=file_now2,into=c("GT","SU","PE","SR"),col=column,sep=":")
		SU_NM<-SU_NM+as.numeric(as.character(column_now$SU))
		SR_NM<-SR_NM+as.numeric(as.character(column_now$SR))
		PE_NM<-PE_NM+as.numeric(as.character(column_now$PE))
		}

	file_now$PE<-test$PE
	file_now$SR<-test$SR
	file_now$SU<-test$SU
	file_now$PE_M<-PE_M
	file_now$SU_M<-SU_M
	file_now$SR_M<-SR_M
	file_now$PE_NM<-PE_NM
	file_now$SU_NM<-SU_NM
	file_now$SR_NM<-SR_NM

	assign(paste0(file,"_DEL"),file_now)

	}

for (file in filelist)
	{file_now1<-get(paste0(file,"_DEL"))
	file_now2<-get(paste0(file,"_DUP"))
	file_now1<-file_now1[((file_now1$PE_NM>10)|(file_now1$PE_M>10)),]
	file_now2<-file_now2[((file_now2$PE_NM>10)|(file_now2$PE_M>10)),]
	assign(paste0(file,"_DEL"),file_now1)
	assign(paste0(file,"_DUP"),file_now2)
	}

require(openxlsx)
for (file in filelist)
	{write.xlsx(get(paste0(file,"_DEL")),paste0(file,"_DEL_Lumpy_Kowa.xlsx"),col.names=TRUE,row.names=FALSE)
	write.xlsx(get(paste0(file,"_DUP")),paste0(file,"_DUP_Lumpy_Kowa.xlsx"),col.names=TRUE,row.names=FALSE)
	}

for (file in filelist)
	{write.table(get(paste0(file,"_DEL")),paste0(file,"_DEL_Lumpy_Kowa.table"),col.names=TRUE,row.names=FALSE,sep="\t")
	write.table(get(paste0(file,"_DUP")),paste0(file,"_DUP_Lumpy_Kowa.table"),col.names=TRUE,row.names=FALSE,sep="\t")
	}

