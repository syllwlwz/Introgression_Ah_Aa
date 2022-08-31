vcf<-read.table("/prj/pflaphy-cutolgs/Filtered_lyrata/BestnewFil2.vcf",sep="\t")
#don't forget to merge in columns later
#Repolarize
repol<-read.table("/prj/pflaphy-cutolgs/Repolarized/Introgression_07repolarized.lookupKey.minInd_4.txt",sep="\t",header=F)
remove<-read.table("/prj/pflaphy-cutolgs/Repolarized/Repol_exclude.txt",sep="\t",header=F)
repol$V1<-as.character(repol$V1)
remove$V1<-as.character(remove$V1)
for (i in c(1,4,5,8:21))
	{vcf[,i]<-as.character(vcf[,i])
	}
repol2<-repol[paste0(repol$V1,repol$V2)%in%paste0(vcf$V1,vcf$V2),]
#90653 out of of 465917

for (i in 1:length(repol))
	{work<-vcf[vcf$V1==repol2[i,1]&vcf$V2==repol2[i,2],]
	Ref<-work$V5
	Alt<-work$V4
	work$V4<-Ref
	work$V5<-Alt
	V8<-strsplit(work$V8,split=";")[[1]]
	AC<-as.numeric(gsub("AC=","",V8[1]))
	AF<-as.numeric(gsub("AF=","",V8[2]))
	AN<-as.numeric(gsub("AN=","",V8[3]))
	AC_new<-AN-AC
	AF_new<-AC_new/AN
	work$V8<-paste0("AC=",AC_new,";AF=",AF_new,";",paste0(V8[3:14],collapse=";"))
	for (j in 10:ncol(work))
		{gtwork<-work[j]
		gtsep<-strsplit(as.character(gtwork),split=":")[[1]]
		gtsep[1]<-gsub("\\|","\\/",gtsep[1])
		gtsep[1]<-ifelse(gtsep[1]=="0/1","1/0",ifelse(gtsep[1]=="0/0","1/1","0/0"))
		gtsep2<-strsplit(gtsep[2],split=",")[[1]]
		gtsep[2]<-paste0(gtsep2[2],",",gtsep2[1])
		gtsep[5]<-ifelse(gtsep[1]=="0|1","1|0",ifelse(gtsep[1]=="0|0","1|1","0|0"))
		gtsep7<-strsplit(gtsep[7],split=",")[[1]]
		gtsep[7]<-paste0(c(gtsep7[3],gtsep7[2],gtsep7[1]),collapse=",")
		work[j]<-paste0(gtsep,collapse=":")
		}
	vcf[vcf$V1==repol2[i,1]&vcf$V2==repol2[i,2],]<-work
	}


vcf2<-vcf[!paste0(vcf$V1,vcf$V2)%in%paste0(remove$V1,remove$V2),]
#2572306 of 2702190
#95%

write.table(vcf2,"/prj/pflaphy-cutolgs/GS_LangBest_fil/Bestnewfil2_repol.vcf",sep="\t",quote=F,row.names=F,col.names=F)

#length(strsplit(gtsep[1],split="/")[[1]]) for polyploids
#need to switch GT,AD,PGT,PL
#GT genotype
#AD allele depth
#GQ genotyping quality
#PGT Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another
#PID "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group"
#PL Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification


