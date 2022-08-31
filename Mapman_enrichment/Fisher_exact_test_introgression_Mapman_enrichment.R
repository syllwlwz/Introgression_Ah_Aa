require(openxlsx)
GS<-read.xlsx("Genome_scans/Selected_genes_Fst_introgression_05_MiasKowa.xlsx",1)
Twisst<-read.xlsx("F:/Introgression/Twisst/Twisst_MiashalleriKowaarenosa_w100/GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_Mias_Twisst_w100_MahalleriKaar_07.xlsx",1)
CNVs_old<-read.xlsx("GROM/CNVs_overlap_strict_MiasKowa_arenosa_Lumpy_dedup_newortho.xlsx",1)
CNVs_cnmopsonly<-read.xlsx("GROM/CNVs_overlap_strict_MiasKowa_arenosa_Lumpy_dedup_newortho_cnmopsonly.xlsx",1)

RNASeq_R<-read.xlsx("Introrna/DESeq2/Deseq2_introgression_inclAth_Ctrl_samedir.xlsx",1)
RNASeq_S<-read.xlsx("Introrna/DESeq2/Deseq2_introgression_inclAth_Ctrl_samedir.xlsx",2)


require(openxlsx)
Mapman<-read.table("F:/Introgression/Mapman_for_merge.table",header=T,sep="\t",fill=T)
Mapman<-unique(Mapman)
Mapman$IDENTIFIER<-toupper(Mapman$IDENTIFIER)
#460 no AGI code but name, with M for metabolite against T for transcript
#str(Mapman[Mapman$TYPE=="M",])
Mapman<-Mapman[Mapman$TYPE=="T",]
require(data.table)
Mapman_cat<-unique(data.frame(Mapman$Main,Mapman$Main_Name))
colnames(Mapman_cat)<-c("Main","Main_Name")
Mapman_cat$Main_Name<-as.character(Mapman_cat$Main_Name)
Mapman_cat$Main_Name[Mapman_cat$Main=="39"]<-"Casparian_stripe_suberin"
Mapman_cat<-Mapman_cat[!duplicated(Mapman_cat[,1:2]),]

OG_Alyrata2<-read.xlsx("F:/Introgression/Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx")
OG_Alyrata3<-OG_Alyrata2[-grep("U",OG_Alyrata2$Alyr_ID),]
OG_Alyrata<-rbind(OG_Alyrata3,OG_Alyrata2[grep("AL9U",OG_Alyrata2$Alyr_ID),])
OG_Alyrata_Twisst<-OG_Alyrata3

OG_Alyrata_MM<-merge(OG_Alyrata,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
OG_Alyrata_MM<-OG_Alyrata_MM[order(OG_Alyrata_MM$Main),]
str(OG_Alyrata_MM)
#26266OG_Alyrata_MM_table2<-table(OG_Alyrata_MM$Main)
OG_Alyrata_MM_table<-merge(as.data.frame(OG_Alyrata_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
OG_Alyrata_MM_table<-OG_Alyrata_MM_table[,c(1,3,2)]
OG_Alyrata_MM_table<-OG_Alyrata_MM_table[order(OG_Alyrata_MM_table$Var1),]
OG_Alyrata_MM_table_n<-nrow(OG_Alyrata_MM)



OG_Alyrata_MM_Twisst<-merge(OG_Alyrata_Twisst,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
OG_Alyrata_MM_Twisst<-OG_Alyrata_MM_Twisst[order(OG_Alyrata_MM_Twisst$Main),]
str(OG_Alyrata_MM_Twisst)
#26147

OG_Alyrata_MM_Twisst_table2<-table(OG_Alyrata_MM_Twisst$Main)
OG_Alyrata_MM_Twisst_table<-merge(as.data.frame(OG_Alyrata_MM_Twisst_table2),Mapman_cat,by.x="Var1",by.y="Main")
OG_Alyrata_MM_Twisst_table<-OG_Alyrata_MM_Twisst_table[,c(1,3,2)]
OG_Alyrata_MM_Twisst_table<-OG_Alyrata_MM_Twisst_table[order(OG_Alyrata_MM_Twisst_table$Var1),]
OG_Alyrata_MM_Twisst_table_n<-nrow(OG_Alyrata_MM_Twisst)


GS_MM<-merge(GS,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
GS_MM<-GS_MM[order(GS_MM$Main),]
GS_MM_table2<-table(GS_MM$Main)
GS_MM_table<-merge(as.data.frame(GS_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
GS_MM_table<-GS_MM_table[,c(1,3,2)]
GS_MM_table<-GS_MM_table[order(GS_MM_table$Var1),]
GS_MM_table_n<-nrow(GS_MM)

x=GS_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=OG_Alyrata_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- GS_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- OG_Alyrata_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(GS_MM_table_n))/(dat[2,1]/OG_Alyrata_MM_table_n))
}

GS_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
GS_fisher$q_value<-qvalue(GS_fisher$p.value,lambda=0)$qvalues
GS_fisher_sig<-GS_fisher[GS_fisher$q_value<0.05,]
write.table(GS_fisher_sig, "GS_fisher.txt",sep="\t",row.names=F,quote=F)


Twisst_MM<-merge(Twisst,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
Twisst_MM<-Twisst_MM[order(Twisst_MM$Main),]
Twisst_MM_table2<-table(Twisst_MM$Main)
Twisst_MM_table<-merge(as.data.frame(Twisst_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
Twisst_MM_table<-Twisst_MM_table[,c(1,3,2)]
Twisst_MM_table<-Twisst_MM_table[order(Twisst_MM_table$Var1),]
Twisst_MM_table_n<-nrow(Twisst_MM)
OG_Alyrata_MM_table_n_Twisst<-nrow(OG_Alyrata_MM_Twisst)

x=Twisst_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=OG_Alyrata_MM_Twisst_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- Twisst_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- OG_Alyrata_MM_Twisst_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(Twisst_MM_table_n))/(dat[2,1]/OG_Alyrata_MM_Twisst_table_n))
}

Twisst_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
Twisst_fisher$q_value<-qvalue(Twisst_fisher$p.value,lambda=0)$qvalues
Twisst_fisher_sig<-Twisst_fisher[Twisst_fisher$q_value<0.05,]
write.table(Twisst_fisher_sig, "Twisst_fisher.txt",sep="\t",row.names=F,quote=F)


CNVs_old_up<-CNVs_old[!(CNVs_old$CN_class=="CN0"|CNVs_old$CN_class=="CN1"|CNVs_old$CN_class=="CN2"|CNVs_old$CN_class=="CN3"),]
CNVs_old_down<-CNVs_old[(CNVs_old$CN_class=="CN0"|CNVs_old$CN_class=="CN1"|CNVs_old$CN_class=="CN2"|CNVs_old$CN_class=="CN3"),]

CNVs_cnmopsonly_up<-CNVs_cnmopsonly[!(CNVs_cnmopsonly$CN_class=="CN0"|CNVs_cnmopsonly$CN_class=="CN1"|CNVs_cnmopsonly$CN_class=="CN2"|CNVs_cnmopsonly$CN_class=="CN3"),]
CNVs_cnmopsonly_down<-CNVs_cnmopsonly[(CNVs_cnmopsonly$CN_class=="CN0"|CNVs_cnmopsonly$CN_class=="CN1"|CNVs_cnmopsonly$CN_class=="CN2"|CNVs_cnmopsonly$CN_class=="CN3"),]

CNVs_old_up_MM<-merge(CNVs_old_up,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
CNVs_old_up_MM<-CNVs_old_up_MM[order(CNVs_old_up_MM$Main),]
CNVs_old_up_MM_table2<-table(CNVs_old_up_MM$Main)
CNVs_old_up_MM_table<-merge(as.data.frame(CNVs_old_up_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
CNVs_old_up_MM_table<-CNVs_old_up_MM_table[,c(1,3,2)]
CNVs_old_up_MM_table<-CNVs_old_up_MM_table[order(CNVs_old_up_MM_table$Var1),]
CNVs_old_up_MM_table_n<-nrow(CNVs_old_up_MM)

x=CNVs_old_up_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=OG_Alyrata_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- CNVs_old_up_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- OG_Alyrata_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(CNVs_old_up_MM_table_n))/(dat[2,1]/OG_Alyrata_MM_table_n))
}

CNVs_old_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
CNVs_old_up_fisher$q_value<-qvalue(CNVs_old_up_fisher$p.value,lambda=0)$qvalues
CNVs_old_up_fisher_sig<-CNVs_old_up_fisher[CNVs_old_up_fisher$q_value<0.05,]
write.table(CNVs_old_up_fisher_sig, "CNVs_old_up_fisher.txt",sep="\t",row.names=F,quote=F)

CNVs_old_down_MM<-merge(CNVs_old_down,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
CNVs_old_down_MM<-CNVs_old_down_MM[order(CNVs_old_down_MM$Main),]
CNVs_old_down_MM_table2<-table(CNVs_old_down_MM$Main)
CNVs_old_down_MM_table<-merge(as.data.frame(CNVs_old_down_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
CNVs_old_down_MM_table<-CNVs_old_down_MM_table[,c(1,3,2)]
CNVs_old_down_MM_table<-CNVs_old_down_MM_table[order(CNVs_old_down_MM_table$Var1),]
CNVs_old_down_MM_table_n<-nrow(CNVs_old_down_MM)

x=CNVs_old_down_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=OG_Alyrata_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- CNVs_old_down_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- OG_Alyrata_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(CNVs_old_down_MM_table_n))/(dat[2,1]/OG_Alyrata_MM_table_n))
}

CNVs_old_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
CNVs_old_down_fisher$q_value<-qvalue(CNVs_old_down_fisher$p.value,lambda=0)$qvalues
CNVs_old_down_fisher_sig<-CNVs_old_down_fisher[CNVs_old_down_fisher$q_value<0.05,]
write.table(CNVs_old_down_fisher_sig, "CNVs_old_down_fisher.txt",sep="\t",row.names=F,quote=F)

CNVs_cnmopsonly_up_MM<-merge(CNVs_cnmopsonly_up,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
CNVs_cnmopsonly_up_MM<-CNVs_cnmopsonly_up_MM[order(CNVs_cnmopsonly_up_MM$Main),]
CNVs_cnmopsonly_up_MM_table2<-table(CNVs_cnmopsonly_up_MM$Main)
CNVs_cnmopsonly_up_MM_table<-merge(as.data.frame(CNVs_cnmopsonly_up_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
CNVs_cnmopsonly_up_MM_table<-CNVs_cnmopsonly_up_MM_table[,c(1,3,2)]
CNVs_cnmopsonly_up_MM_table<-CNVs_cnmopsonly_up_MM_table[order(CNVs_cnmopsonly_up_MM_table$Var1),]
CNVs_cnmopsonly_up_MM_table_n<-nrow(CNVs_cnmopsonly_up_MM)

x=CNVs_cnmopsonly_up_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=OG_Alyrata_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- CNVs_cnmopsonly_up_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- OG_Alyrata_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(CNVs_cnmopsonly_up_MM_table_n))/(dat[2,1]/OG_Alyrata_MM_table_n))
}

CNVs_cnmopsonly_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
CNVs_cnmopsonly_up_fisher$q_value<-qvalue(CNVs_cnmopsonly_up_fisher$p.value,lambda=0)$qvalues
CNVs_cnmopsonly_up_fisher_sig<-CNVs_cnmopsonly_up_fisher[CNVs_cnmopsonly_up_fisher$q_value<0.05,]
write.table(CNVs_cnmopsonly_up_fisher_sig, "CNVs_cnmopsonly_up_fisher.txt",sep="\t",row.names=F,quote=F)


CNVs_cnmopsonly_down_MM<-merge(CNVs_cnmopsonly_down,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
CNVs_cnmopsonly_down_MM<-CNVs_cnmopsonly_down_MM[order(CNVs_cnmopsonly_down_MM$Main),]
CNVs_cnmopsonly_down_MM_table2<-table(CNVs_cnmopsonly_down_MM$Main)
CNVs_cnmopsonly_down_MM_table<-merge(as.data.frame(CNVs_cnmopsonly_down_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
CNVs_cnmopsonly_down_MM_table<-CNVs_cnmopsonly_down_MM_table[,c(1,3,2)]
CNVs_cnmopsonly_down_MM_table<-CNVs_cnmopsonly_down_MM_table[order(CNVs_cnmopsonly_down_MM_table$Var1),]
CNVs_cnmopsonly_down_MM_table_n<-nrow(CNVs_cnmopsonly_down_MM)

x=CNVs_cnmopsonly_down_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=OG_Alyrata_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- CNVs_cnmopsonly_down_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- OG_Alyrata_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(CNVs_cnmopsonly_down_MM_table_n))/(dat[2,1]/OG_Alyrata_MM_table_n))
}

CNVs_cnmopsonly_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
CNVs_cnmopsonly_down_fisher$q_value<-qvalue(CNVs_cnmopsonly_down_fisher$p.value,lambda=0)$qvalues
CNVs_cnmopsonly_down_fisher_sig<-CNVs_cnmopsonly_down_fisher[CNVs_cnmopsonly_down_fisher$q_value<0.05,]
write.table(CNVs_cnmopsonly_down_fisher_sig, "CNVs_cnmopsonly_down_fisher.txt",sep="\t",row.names=F,quote=F)


coldata<-read.csv("Introrna/DESeq2/Introrna_sample_info.csv",sep=";")
FPKM<-read.table("Introrna/DESeq2/FPKM_df.txt",sep="\t",header=T,row.names=1)
CPM_df<-FPKM
for (i in 1:60)
	{Si<-data.frame(CPM_df[,1+i]/(sum(CPM_df[,1+i])/1000000))
	colnames(Si)<-paste("CPM_S",i,sep="")
	CPM_df<-data.frame(CPM_df,Si)
	}
Log_cpm<-log2(CPM_df[,122:181]+1)
colnames(Log_cpm)<-apply(coldata[,2:6],1,function(x) paste(x,sep="",collapse="_"))
rownames(Log_cpm)<-rownames(FPKM)
Log_cpm_Ctrl<-Log_cpm[,grep("c_",colnames(Log_cpm))]


#root
All<-merge(OG_Alyrata2,Log_cpm_Ctrl,by.y="row.names",by.x="Alyr_ID")
All_cpm_R<-All[apply(All[,grep("_R_",colnames(All))],1,function(x) any(x>1.4)),]

All_R_MM<-merge(All_cpm_R,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
All_R_MM<-All_R_MM[order(All_R_MM$Main),]
str(All_R_MM)
#21986
#dupl Ath IDs because some in several categories
#All_R_MM<-All_R_MM[!duplicated(All_R_MM$Lyr_Gene),]

All_R_MM_table2<-table(All_R_MM$Main)
All_R_MM_table<-merge(as.data.frame(All_R_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
All_R_MM_table<-All_R_MM_table[,c(1,3,2)]
All_R_MM_table<-All_R_MM_table[order(All_R_MM_table$Var1),]
All_R_MM_table_n<-nrow(All_R_MM)

#shoot
All<-merge(OG_Alyrata2,Log_cpm_Ctrl,by.y="row.names",by.x="Alyr_ID")
All_cpm_S<-All[apply(All[,grep("_S_",colnames(All))],1,function(x) any(x>1.4)),]

All_S_MM<-merge(All_cpm_S,Mapman,by.x="Ath_ID",by.y="IDENTIFIER")
All_S_MM<-All_S_MM[order(All_S_MM$Main),]
str(All_S_MM)
#21936
#dupl Ath IDs because some in several categories
#All_S_MM<-All_S_MM[!duplicated(All_S_MM$Ah_ID),]

All_S_MM_table2<-table(All_S_MM$Main)
All_S_MM_table<-merge(as.data.frame(All_S_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
All_S_MM_table<-All_S_MM_table[,c(1,3,2)]
All_S_MM_table<-All_S_MM_table[order(All_S_MM_table$Var1),]
All_S_MM_table_n<-nrow(All_S_MM)

RNASeq_S_cpm<-merge(RNASeq_S,Log_cpm,by.y="row.names",by.x="Alyr_ID")
RNASeq_S_cpm_fil<-RNASeq_S_cpm[apply(RNASeq_S_cpm[,grep("c_S_",colnames(RNASeq_S_cpm[,326:385]))+325],1,function(x) any(x>1.4)),]

RNASeq_R_cpm<-merge(RNASeq_R,Log_cpm,by.y="row.names",by.x="Alyr_ID")
RNASeq_R_cpm_fil<-RNASeq_R_cpm[apply(RNASeq_R_cpm[,grep("c_R_",colnames(RNASeq_S_cpm[,326:385]))+325],1,function(x) any(x>1.4)),]

RNASeq_R_cpm_Mapman_up<-RNASeq_R_cpm_fil[RNASeq_R_cpm_fil$ log2FoldChange.x>0,]
RNASeq_R_cpm_Mapman_down<-RNASeq_R_cpm_fil[RNASeq_R_cpm_fil$ log2FoldChange.x<0,]

RNASeq_S_cpm_Mapman_up<-RNASeq_S_cpm_fil[RNASeq_S_cpm_fil$ log2FoldChange.x>0,]
RNASeq_S_cpm_Mapman_down<-RNASeq_S_cpm_fil[RNASeq_S_cpm_fil$ log2FoldChange.x<0,]

RNASeq_S_up_MM<-merge(RNASeq_S_cpm_Mapman_up,Mapman,by.x="Ath_ID.x",by.y="IDENTIFIER")
RNASeq_S_up_MM<-RNASeq_S_up_MM[order(RNASeq_S_up_MM$Main),]
RNASeq_S_up_MM_table2<-table(RNASeq_S_up_MM$Main)
RNASeq_S_up_MM_table<-merge(as.data.frame(RNASeq_S_up_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
RNASeq_S_up_MM_table<-RNASeq_S_up_MM_table[,c(1,3,2)]
RNASeq_S_up_MM_table<-RNASeq_S_up_MM_table[order(RNASeq_S_up_MM_table$Var1),]
RNASeq_S_up_MM_table_n<-nrow(RNASeq_S_up_MM)

x=RNASeq_S_up_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=All_S_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- RNASeq_S_up_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- All_S_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(RNASeq_S_up_MM_table_n))/(dat[2,1]/All_S_MM_table_n))
}

RNASeq_S_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
RNASeq_S_up_fisher$q_value<-qvalue(RNASeq_S_up_fisher$p.value,lambda=0)$qvalues
RNASeq_S_up_fisher_sig<-RNASeq_S_up_fisher[RNASeq_S_up_fisher$q_value<0.05,]
write.table(RNASeq_S_up_fisher_sig, "RNASeq_S_up_fisher.txt",sep="\t",row.names=F,quote=F)


RNASeq_S_down_MM<-merge(RNASeq_S_cpm_Mapman_down,Mapman,by.x="Ath_ID.x",by.y="IDENTIFIER")
RNASeq_S_down_MM<-RNASeq_S_down_MM[order(RNASeq_S_down_MM$Main),]
RNASeq_S_down_MM_table2<-table(RNASeq_S_down_MM$Main)
RNASeq_S_down_MM_table<-merge(as.data.frame(RNASeq_S_down_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
RNASeq_S_down_MM_table<-RNASeq_S_down_MM_table[,c(1,3,2)]
RNASeq_S_down_MM_table<-RNASeq_S_down_MM_table[order(RNASeq_S_down_MM_table$Var1),]
RNASeq_S_down_MM_table_n<-nrow(RNASeq_S_down_MM)

x=RNASeq_S_down_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=All_S_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- RNASeq_S_down_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- All_S_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(RNASeq_S_down_MM_table_n))/(dat[2,1]/All_S_MM_table_n))
}

RNASeq_S_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
RNASeq_S_down_fisher$q_value<-qvalue(RNASeq_S_down_fisher$p.value,lambda=0)$qvalues
RNASeq_S_down_fisher_sig<-RNASeq_S_down_fisher[RNASeq_S_down_fisher$q_value<0.05,]
#write.table(RNASeq_S_down_fisher_sig, "RNASeq_S_down_fisher.txt",sep="\t",row.names=F,quote=F)



RNASeq_R_up_MM<-merge(RNASeq_R_cpm_Mapman_up,Mapman,by.x="Ath_ID.x",by.y="IDENTIFIER")
RNASeq_R_up_MM<-RNASeq_R_up_MM[order(RNASeq_R_up_MM$Main),]
RNASeq_R_up_MM_table2<-table(RNASeq_R_up_MM$Main)
RNASeq_R_up_MM_table<-merge(as.data.frame(RNASeq_R_up_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
RNASeq_R_up_MM_table<-RNASeq_R_up_MM_table[,c(1,3,2)]
RNASeq_R_up_MM_table<-RNASeq_R_up_MM_table[order(RNASeq_R_up_MM_table$Var1),]
RNASeq_R_up_MM_table_n<-nrow(RNASeq_R_up_MM)

x=RNASeq_R_up_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=All_R_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- RNASeq_R_up_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- All_R_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(RNASeq_R_up_MM_table_n))/(dat[2,1]/All_R_MM_table_n))
}

RNASeq_R_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
RNASeq_R_up_fisher$q_value<-qvalue(RNASeq_R_up_fisher$p.value,lambda=0)$qvalues
RNASeq_R_up_fisher_sig<-RNASeq_R_up_fisher[RNASeq_R_up_fisher$q_value<0.05,]
write.table(RNASeq_R_up_fisher_sig, "RNASeq_R_up_fisher.txt",sep="\t",row.names=F,quote=F)


RNASeq_R_down_MM<-merge(RNASeq_R_cpm_Mapman_down,Mapman,by.x="Ath_ID.x",by.y="IDENTIFIER")
RNASeq_R_down_MM<-RNASeq_R_down_MM[order(RNASeq_R_down_MM$Main),]
RNASeq_R_down_MM_table2<-table(RNASeq_R_down_MM$Main)
RNASeq_R_down_MM_table<-merge(as.data.frame(RNASeq_R_down_MM_table2),Mapman_cat,by.x="Var1",by.y="Main")
RNASeq_R_down_MM_table<-RNASeq_R_down_MM_table[,c(1,3,2)]
RNASeq_R_down_MM_table<-RNASeq_R_down_MM_table[order(RNASeq_R_down_MM_table$Var1),]
RNASeq_R_down_MM_table_n<-nrow(RNASeq_R_down_MM)

x=RNASeq_R_down_MM_table[1:3]
x[,3]<-as.numeric(as.character(x[,3]))

y=All_R_MM_table
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$Var1)) 
{  
  dat[1,1] <- x$Freq[x$Var1==i]
  dat[1,2] <- RNASeq_R_down_MM_table_n-dat[1,1]
  dat[2,1] <- y$Freq[y$Var1==i]-dat[1,1]
  dat[2,2] <- All_R_MM_table_n-y$Freq[y$Var1==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(RNASeq_R_down_MM_table_n))/(dat[2,1]/All_R_MM_table_n))
}

RNASeq_R_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
RNASeq_R_down_fisher$q_value<-qvalue(RNASeq_R_down_fisher$p.value,lambda=0)$qvalues
RNASeq_R_down_fisher_sig<-RNASeq_R_down_fisher[RNASeq_R_down_fisher$q_value<0.05,]
write.table(RNASeq_R_down_fisher_sig, "RNASeq_R_down_fisher.txt",sep="\t",row.names=F,quote=F)
#positive log2fold change = Zapa higher expression



bpRNASeq_R_down_fisher_sig<-barplot(sort(-log10(RNASeq_R_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(RNASeq_R_down_fisher_sig$q_value))),width=1.2)
bpRNASeq_R_up_fisher_sig<-barplot(sort(-log10(RNASeq_R_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(RNASeq_R_up_fisher_sig$q_value))),width=1.2)
#bpRNASeq_S_down_fisher_sig<-barplot(sort(-log10(RNASeq_S_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(RNASeq_S_down_fisher_sig$q_value))),width=1.2)
bpRNASeq_S_up_fisher_sig<-barplot(sort(-log10(RNASeq_S_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(RNASeq_S_up_fisher_sig$q_value))),width=1.2)

pdf('GOenrichment_Mapman_RNASeq_introgression_Ctrl_upanddownsep_qvalue.pdf', width=9, height=6,paper="special",pointsize=15)
par(oma=c(6,4,2,1))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),heights=c(1,0.75),widths=c(1,0.84))
par(mar=c(2,10,1,0),xpd=T)

barplot(sort(-log10(RNASeq_R_up_fisher_sig$q_value)),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,15),col="blue",cex.axis=1.3,axes=F)
axis(side=2,tick=F,labels=c("Metal\nhomeostasis","Secondary\nmetabolism","Misc"),at=bpRNASeq_R_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
#mtext(text=expression("Transcript abundance higher in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.3)

barplot(sort(-log10(RNASeq_S_up_fisher_sig$q_value)),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,5),col="blue",cex.axis=1.3,axes=F)
axis(side=2,tick=F,labels=c("Stress","Secondary\nmetabolism","Misc."),at=bpRNASeq_S_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(-log[10](adjusted~italic(p)-value)),side=1,cex=1.3,outer=T,line=1)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(-log10(RNASeq_R_down_fisher_sig$q_value)),beside=T,horiz=T,las=1,cex.lab=1.5,xlim=c(0,15),width=1.2,col="red",cex.axis=1.3)
axis(side=2,tick=F,labels=c("Metal\nhomeostasis","Stress"),at=bpRNASeq_R_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(Enriched~Mapman~categories),side=2,cex=1.3,outer=T,line=2,xlim=c(0,15))
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(0,0,xlim=c(0,5),cex.axis=1.3,beside=T,las=1,horiz=T)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), cnmopsonly = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend=c(expression(paste("Transcript abundance lower in Mias ",italic(A.~arenosa)," than Zapa"~italic(A.~arenosa)),paste("Transcript abundance higher in Mias ",italic(A.~arenosa)," than Zapa ",italic(A.~arenosa)))),fill=c("blue","red"),xpd=T,bty="n",inset=c(0,0),cex=1.3)


dev.off()


bpRNASeq_R_down_fisher_sig<-barplot(sort(-log10(RNASeq_R_down_fisher_sig$Fold_enrichment)),las=3,xlim=c(0,max(-log10(RNASeq_R_down_fisher_sig$Fold_enrichment))),width=1.2)
bpRNASeq_R_up_fisher_sig<-barplot(sort(-log10(RNASeq_R_up_fisher_sig$Fold_enrichment)),las=3,xlim=c(0,max(-log10(RNASeq_R_up_fisher_sig$Fold_enrichment))),width=1.2)
#bpRNASeq_S_down_fisher_sig<-barplot(sort(-log10(RNASeq_S_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(RNASeq_S_down_fisher_sig$q_value))),width=1.2)
bpRNASeq_S_up_fisher_sig<-barplot(sort(-log10(RNASeq_S_up_fisher_sig$Fold_enrichment)),las=3,xlim=c(0,max(-log10(RNASeq_S_up_fisher_sig$Fold_enrichment))),width=1.2)


pdf('GOenrichment_Mapman_RNASeq_introgression_Ctrl_upanddownsep_foldenrichment.pdf', width=9, height=6,paper="special",pointsize=15)
par(oma=c(6,4,2,1))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),heights=c(1,0.75),widths=c(1,0.84))
par(mar=c(2,10,1,0),xpd=T)

barplot(sort(RNASeq_R_up_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,5),col="blue",cex.axis=1.3,axes=F)
axis(side=2,tick=F,labels=c("Misc.","Secondary\nmetabolism","Metal\nhomeostasis"),at=bpRNASeq_R_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
#mtext(text=expression("Transcript abundance higher in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.5)

barplot(sort(RNASeq_S_up_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,5),col="blue",cex.axis=1.3,axes=F)
axis(side=2,tick=F,labels=c("Stress","Misc.","Secondary\nmetabolism"),at=bpRNASeq_S_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text="Fold enrichment",side=1,cex=1.3,outer=T,line=1)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(RNASeq_R_down_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,xlim=c(0,5),width=1.2,col="red",cex.axis=1.3)
axis(side=2,tick=F,labels=c("Stress","Metal\nhomeostasis"),at=bpRNASeq_R_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(Enriched~Mapman~categories),side=2,cex=1.3,outer=T,line=2,xlim=c(0,15))
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(0,0,xlim=c(0,5),cex.axis=1.3,beside=T,las=1,horiz=T)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), cnmopsonly = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend=c(expression(paste("Transcript abundance lower in Mias ",italic(A.~arenosa)," than Zapa"~italic(A.~arenosa)),paste("Transcript abundance higher in Mias ",italic(A.~arenosa)," than Zapa ",italic(A.~arenosa)))),fill=c("blue","red"),xpd=T,bty="n",inset=c(0,0),cex=1.3)


dev.off()













Twisst_fisher_sig
CNVs_old_down_fisher_sig
CNVs_old_up_fisher_sig
CNVs_cnmopsonly_down_fisher_sig
CNVs_cnmopsonly_up_fisher_sig
GS_fisher_sig
RNASeq_R_down_fisher_sig
RNASeq_R_up_fisher_sig
RNASeq_S_down_fisher_sig
RNASeq_S_up_fisher_sig



bpTwisst_fisher_sig<-barplot(sort(Twisst_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(Twisst_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
#bpCNVs_old_down_fisher_sig<-barplot(sort(CNVs_old_down_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(CNVs_old_down_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
bpCNVs_old_up_fisher_sig<-barplot(sort(CNVs_old_up_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(CNVs_old_up_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
bpCNVs_cnmopsonly_down_fisher_sig<-barplot(sort(CNVs_cnmopsonly_down_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(CNVs_cnmopsonly_down_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
bpCNVs_cnmopsonly_up_fisher_sig<-barplot(sort(CNVs_cnmopsonly_up_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(CNVs_cnmopsonly_up_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
bpGS_fisher_sig<-barplot(sort(GS_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(GS_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
bpRNASeq_R_down_fisher_sig<-barplot(sort(RNASeq_R_down_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(RNASeq_R_down_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
bpRNASeq_R_up_fisher_sig<-barplot(sort(RNASeq_R_up_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(RNASeq_R_up_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
#bpRNASeq_S_down_fisher_sig<-barplot(sort(RNASeq_S_down_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(RNASeq_S_down_fisher_sig$Fold_enrichment)),width=1.2,plot=F)
bpRNASeq_S_up_fisher_sig<-barplot(sort(RNASeq_S_up_fisher_sig$Fold_enrichment),las=3,xlim=c(0,max(RNASeq_S_up_fisher_sig$Fold_enrichment)),width=1.2,plot=F)

Twisst_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",Twisst_fisher_sig$Main_Name))))))))))
CNVs_old_up_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",CNVs_old_up_fisher_sig$Main_Name))))))))))
CNVs_cnmopsonly_down_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",CNVs_cnmopsonly_down_fisher_sig$Main_Name))))))))))
CNVs_cnmopsonly_up_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",CNVs_cnmopsonly_up_fisher_sig$Main_Name))))))))))
GS_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",GS_fisher_sig$Main_Name))))))))))
RNASeq_R_down_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",RNASeq_R_down_fisher_sig$Main_Name))))))))))
RNASeq_R_up_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",RNASeq_R_up_fisher_sig$Main_Name))))))))))
RNASeq_S_up_fisher_sig$Main_Name<-gsub("secondary metabolism","Secondary metabolism",gsub("misc","Misc.",gsub("PS","Photosynthesis",gsub("protein","Protein",gsub("FE","Fe-sulfur cluster",gsub("transport","Transport",gsub("metal","Metal\nhomeostasis",gsub("stress","Stress",gsub("nutrition&ionome","Nutrition&ionome",gsub("mitochondrial electron transport / ATP synthesis","Mitochondrial electron transport/\n ATP synthesis",RNASeq_S_up_fisher_sig$Main_Name))))))))))



pdf('GOenrichment_Mapman_introgression_all_Twisst_Fold_enrichment.pdf',width=18,height=8,paper="special",pointsize=22)

par(oma=c(1,3,1,3))
layout(mat=matrix(c(1,2,3,4,5,6),nrow=3,ncol=2),heights=c(0.8,0.8,0.8),widths=c(1,0.75))
par(mar=c(1,15,2,1),xpd=T)

barplot(sort(Twisst_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,20),col=rgb(0.2,1,0.2,alpha=0.3),cex.axis=1.5,main="Introgression",cex.main=1.5,names=Twisst_fisher_sig$Main_Name,cex.names=1.3,axes=F)
barplot(sort(CNVs_old_up_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,20),col=rgb(0.1,0.4,1,alpha=0.3),cex.axis=1.5,main="Copy number expansion",cex.main=1.5,names=CNVs_old_up_fisher_sig$Main_Name,cex.names=1.3,axes=F)
barplot(sort(GS_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,20),col=rgb(1,0.2,0.2,alpha=0.3),cex.axis=1.5,main="Selection",cex.main=1.5,names=GS_fisher_sig$Main_Name,cex.names=1.3)

barplot(sort(RNASeq_R_down_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,5),col=rgb(1,0.6,0,alpha=0.3),cex.axis=1.5,main="Root Mias higher expression",cex.main=1.5,names=RNASeq_R_down_fisher_sig$Main_Name,cex.names=1.3,axes=F)
barplot(sort(RNASeq_R_up_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,5),col=rgb(1,0.6,0,alpha=0.3),cex.axis=1.5,main="Root Mias lower expression",cex.main=1.5,names=RNASeq_R_up_fisher_sig$Main_Name,cex.names=1.3,axes=F)
barplot(sort(RNASeq_S_up_fisher_sig$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,5),col=rgb(1,0.6,0,alpha=0.3),cex.axis=1.5,main="Shoot Mias lower expression",cex.main=1.5,names=RNASeq_S_up_fisher_sig$Main_Name,cex.names=1.3)

dev.off()


