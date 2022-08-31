##########################################################
#Load and format coverage data per gene for lyrata genome#
##########################################################
options(java.parameters = "-Xmx20000m")

library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)
require(openxlsx)

Thal_description<-read.xlsx("F:/Lyrata_RBH/Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx")

UtesMetal=read.xlsx("../metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

###########################
#Kowa as outgroup#
###########################
############################################################
#GROM data#
############################################################
Buko_arenosa_CNVs1<-read.table("Buko.arenosa.cnvs.vcf",sep="\t",header=F)

names(Buko_arenosa_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Buko_arenosa_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Buko_arenosa_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Buko_arenosa_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Buko_arenosa_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Buko_arenosa_CNVs<-data.frame(Buko_arenosa_CNVs1[,1:2],as.numeric(as.character(Buko_arenosa_CNVs_END[,2])),Buko_arenosa_CNVs1[,3:7],Buko_arenosa_CNVs_info)
names(Buko_arenosa_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Kowa_arenosa_CNVs1<-read.table("Kowa.arenosa.cnvs.vcf",sep="\t",header=F)

names(Kowa_arenosa_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Kowa_arenosa_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_arenosa_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Kowa_arenosa_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_arenosa_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Kowa_arenosa_CNVs<-data.frame(Kowa_arenosa_CNVs1[,1:2],as.numeric(as.character(Kowa_arenosa_CNVs_END[,2])),Kowa_arenosa_CNVs1[,3:7],Kowa_arenosa_CNVs_info)
names(Kowa_arenosa_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Buko_arenosa_cnvs_HC_strict<-read.table("../Cnmops/Kowa/BKar_cnvs_only_strict_minL6_WL350.table",header=T)
Kowa_arenosa_cnvs_HC_strict<-read.table("../Cnmops/Kowa/KBar_cnvs_only_strict_minL6_WL350.table",header=T)

Buko_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Buko_arenosa_CNVs$Scaffold),ranges=IRanges(start=Buko_arenosa_CNVs$Start_pos,end=Buko_arenosa_CNVs$End_pos))
values(Buko_arenosa_CNVs_GRange)<-Buko_arenosa_CNVs[,4:12]
Kowa_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_arenosa_CNVs$Scaffold),ranges=IRanges(start=Kowa_arenosa_CNVs$Start_pos,end=Kowa_arenosa_CNVs$End_pos))
values(Kowa_arenosa_CNVs_GRange)<-Kowa_arenosa_CNVs[,4:12]
BKar_merge<-mergeByOverlaps(Buko_arenosa_CNVs_GRange,Kowa_arenosa_CNVs_GRange)
BKar_merge_same<-BKar_merge[abs(as.numeric(as.character(BKar_merge@ listData$ Buko_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(BKar_merge@ listData$ Kowa_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
BKar_merge_same_df<-as.data.frame(BKar_merge_same)[,1:3]
names(BKar_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
BKar_merge_same_df_GRange<-GRanges(seqnames=tolower(BKar_merge_same_df$Scaffold),ranges=IRanges(start=BKar_merge_same_df$Start_pos,end=BKar_merge_same_df$End_pos))

#subset Buko_arenosa CNVs to regions without similar CN in Kowa_arenosa

BKar_diff<-setdiff(Buko_arenosa_CNVs_GRange,BKar_merge_same_df_GRange)
BKar<-subsetByOverlaps(Buko_arenosa_CNVs_GRange,BKar_diff)

#overlap cn.mops and GROM regions

Buko_arenosa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Buko_arenosa_cnvs_HC_strict$Chr),ranges=IRanges(start=Buko_arenosa_cnvs_HC_strict$gene_start,end=Buko_arenosa_cnvs_HC_strict$gene_end))
values(Buko_arenosa_cnvs_HC_strict_GRange)<-cbind(Buko_arenosa_cnvs_HC_strict[,1:12],Buko_arenosa_cnvs_HC_strict[,18:40])
Buko_arenosa_overlap_strict<-mergeByOverlaps(Buko_arenosa_cnvs_HC_strict_GRange,BKar)
Buko_arenosa_overlap_strict_df<-as.data.frame(Buko_arenosa_overlap_strict)

Buko_arenosa_CNVs_overlapping_strict<-data.frame(Buko_arenosa_overlap_strict_df[,1:40],Buko_arenosa_overlap_strict_df[,76:89])
names(Buko_arenosa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Buko_arenosa_CNVs_HC_GROM_strict<-Buko_arenosa_CNVs_overlapping_strict

CNVs_Buko_arenosa_HC_GROM_cov<-Buko_arenosa_CNVs_HC_GROM_strict
CNVs_Buko_arenosa_HC_GROM_cov<-CNVs_Buko_arenosa_HC_GROM_cov[!((CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN3")&CNVs_Buko_arenosa_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Buko_arenosa_HC_GROM_cov<-CNVs_Buko_arenosa_HC_GROM_cov[!(CNVs_Buko_arenosa_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Buko_arenosa_HC_GROM_cov$CN_class=="CN3")),]

write.table(CNVs_Buko_arenosa_HC_GROM_cov,"Buko_Kowa_arenosa_CNVs_overlap_strict.table",row.names=F,sep="\t")
require(openxlsx)
write.xlsx(CNVs_Buko_arenosa_HC_GROM_cov,"CNVs_overlap_strict_BukoKowa_arenosa.xlsx",sheetName="Buko_Kowa_arenosa",row.names=F,overwrite=T)

#merge in Lumpy data
require(tidyr)
BZ_ar_Lumpy_DUP<-read.xlsx("../Lumpy/Buko_Kowa_arenosa.vcf_DUP_Lumpy.xlsx",1)
BZ_ar_Lumpy_DEL<-read.xlsx("../Lumpy/Buko_Kowa_arenosa.vcf_DEL_Lumpy.xlsx",1)
test1<-separate(data=BZ_ar_Lumpy_DUP,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
BZ_ar_Lumpy_DUP$End<-as.numeric(gsub("END=","",test1$NA4))
test2<-separate(data=BZ_ar_Lumpy_DEL,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
BZ_ar_Lumpy_DEL$End<-as.numeric(gsub("END=","",test2$NA4))
BZ_ar_Lumpy<-rbind(BZ_ar_Lumpy_DUP,BZ_ar_Lumpy_DEL)
BZ_ar_Lumpy<-BZ_ar_Lumpy[(BZ_ar_Lumpy$End-BZ_ar_Lumpy$Pos)<100000,]
BZ_ar_Lumpy$PE_M_std<-BZ_ar_Lumpy$PE_M/(BZ_ar_Lumpy$End-BZ_ar_Lumpy$Pos)
BZ_ar_Lumpy$SU_M_std<-BZ_ar_Lumpy$SU_M/(BZ_ar_Lumpy$End-BZ_ar_Lumpy$Pos)
BZ_ar_Lumpy$SR_M_std<-BZ_ar_Lumpy$SR_M/(BZ_ar_Lumpy$End-BZ_ar_Lumpy$Pos)
BZ_ar_Lumpy$PE_NM_std<-BZ_ar_Lumpy$PE_NM/(BZ_ar_Lumpy$End-BZ_ar_Lumpy$Pos)
BZ_ar_Lumpy$SU_NM_std<-BZ_ar_Lumpy$SU_NM/(BZ_ar_Lumpy$End-BZ_ar_Lumpy$Pos)
BZ_ar_Lumpy$SR_NM_std<-BZ_ar_Lumpy$SR_NM/(BZ_ar_Lumpy$End-BZ_ar_Lumpy$Pos)

BZ_ar_Lumpy_GRange<-GRanges(seqnames=tolower(BZ_ar_Lumpy$Scaffold),ranges=IRanges(start=BZ_ar_Lumpy$Pos,end=BZ_ar_Lumpy$End))
values(BZ_ar_Lumpy_GRange)<-cbind(BZ_ar_Lumpy[,8:ncol(BZ_ar_Lumpy)])

CNVs_Buko_arenosa_HC_GROM_cov_GRange<-GRanges(seqnames=tolower(CNVs_Buko_arenosa_HC_GROM_cov$Scaffold),ranges=IRanges(start=ifelse(CNVs_Buko_arenosa_HC_GROM_cov$Copy_start<CNVs_Buko_arenosa_HC_GROM_cov$GROM_Start_pos,CNVs_Buko_arenosa_HC_GROM_cov$Copy_start,CNVs_Buko_arenosa_HC_GROM_cov$GROM_Start_pos),end=ifelse(CNVs_Buko_arenosa_HC_GROM_cov$Copy_end>CNVs_Buko_arenosa_HC_GROM_cov$GROM_End_pos,CNVs_Buko_arenosa_HC_GROM_cov$Copy_end,CNVs_Buko_arenosa_HC_GROM_cov$GROM_End_pos)))
values(CNVs_Buko_arenosa_HC_GROM_cov_GRange)<-cbind(CNVs_Buko_arenosa_HC_GROM_cov[,1:54])

Buko_arenosa_overlap_strict2<-mergeByOverlaps(CNVs_Buko_arenosa_HC_GROM_cov_GRange,BZ_ar_Lumpy_GRange)
Buko_arenosa_overlap_strict_df2<-as.data.frame(Buko_arenosa_overlap_strict2)[,c(60:113,(113+ncol(BZ_ar_Lumpy)-1):((113+ncol(BZ_ar_Lumpy)-1)+(ncol(BZ_ar_Lumpy)-8)))]

write.xlsx(Buko_arenosa_overlap_strict_df2,"CNVs_overlap_strict_BukoKowa_arenosa_Lumpy.xlsx",sheetName="Buko_Kowa_arenosa",row.names=F,overwrite=T)

############################################################
#GROM data#
############################################################
Buko_halleri_CNVs1<-read.table("Buko.halleri.cnvs.vcf",sep="\t",header=F)

names(Buko_halleri_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Buko_halleri_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Buko_halleri_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Buko_halleri_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Buko_halleri_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Buko_halleri_CNVs<-data.frame(Buko_halleri_CNVs1[,1:2],as.numeric(as.character(Buko_halleri_CNVs_END[,2])),Buko_halleri_CNVs1[,3:7],Buko_halleri_CNVs_info)
names(Buko_halleri_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Kowa_halleri_CNVs1<-read.table("Kowa.halleri.cnvs.vcf",sep="\t",header=F)

names(Kowa_halleri_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Kowa_halleri_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_halleri_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Kowa_halleri_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_halleri_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Kowa_halleri_CNVs<-data.frame(Kowa_halleri_CNVs1[,1:2],as.numeric(as.character(Kowa_halleri_CNVs_END[,2])),Kowa_halleri_CNVs1[,3:7],Kowa_halleri_CNVs_info)
names(Kowa_halleri_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Buko_halleri_cnvs_HC_strict<-read.table("../Cnmops/Kowa/BZha_cnvs_only_strict_minL6_WL350.table",header=T)
Kowa_halleri_cnvs_HC_strict<-read.table("../Cnmops/Kowa/ZBha_cnvs_only_strict_minL6_WL350.table",header=T)

Buko_halleri_CNVs_GRange<-GRanges(seqnames=tolower(Buko_halleri_CNVs$Scaffold),ranges=IRanges(start=Buko_halleri_CNVs$Start_pos,end=Buko_halleri_CNVs$End_pos))
values(Buko_halleri_CNVs_GRange)<-Buko_halleri_CNVs[,4:12]
Kowa_halleri_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_halleri_CNVs$Scaffold),ranges=IRanges(start=Kowa_halleri_CNVs$Start_pos,end=Kowa_halleri_CNVs$End_pos))
values(Kowa_halleri_CNVs_GRange)<-Kowa_halleri_CNVs[,4:12]
BZha_merge<-mergeByOverlaps(Buko_halleri_CNVs_GRange,Kowa_halleri_CNVs_GRange)
BZha_merge_same<-BZha_merge[abs(as.numeric(as.character(BZha_merge@ listData$ Buko_halleri_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(BZha_merge@ listData$ Kowa_halleri_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
BZha_merge_same_df<-as.data.frame(BZha_merge_same)[,1:3]
names(BZha_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
BZha_merge_same_df_GRange<-GRanges(seqnames=tolower(BZha_merge_same_df$Scaffold),ranges=IRanges(start=BZha_merge_same_df$Start_pos,end=BZha_merge_same_df$End_pos))

#subset Buko_halleri CNVs to regions without similar CN in Kowa_halleri

BZha_diff<-setdiff(Buko_halleri_CNVs_GRange,BZha_merge_same_df_GRange)
BZha<-subsetByOverlaps(Buko_halleri_CNVs_GRange,BZha_diff)

#overlap cn.mops and GROM regions

Buko_halleri_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Buko_halleri_cnvs_HC_strict$Chr),ranges=IRanges(start=Buko_halleri_cnvs_HC_strict$gene_start,end=Buko_halleri_cnvs_HC_strict$gene_end))
values(Buko_halleri_cnvs_HC_strict_GRange)<-cbind(Buko_halleri_cnvs_HC_strict[,1:12],Buko_halleri_cnvs_HC_strict[,18:40])
Buko_halleri_overlap_strict<-mergeByOverlaps(Buko_halleri_cnvs_HC_strict_GRange,BZha)
Buko_halleri_overlap_strict_df<-as.data.frame(Buko_halleri_overlap_strict)

Buko_halleri_CNVs_overlapping_strict<-data.frame(Buko_halleri_overlap_strict_df[,1:40],Buko_halleri_overlap_strict_df[,76:89])
names(Buko_halleri_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Buko_halleri_CNVs_HC_GROM_strict<-Buko_halleri_CNVs_overlapping_strict

CNVs_Buko_halleri_HC_GROM_cov<-Buko_halleri_CNVs_HC_GROM_strict
CNVs_Buko_halleri_HC_GROM_cov<-CNVs_Buko_halleri_HC_GROM_cov[!((CNVs_Buko_halleri_HC_GROM_cov$CN_class=="CN0"|CNVs_Buko_halleri_HC_GROM_cov$CN_class=="CN1")&CNVs_Buko_halleri_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Buko_halleri_HC_GROM_cov<-CNVs_Buko_halleri_HC_GROM_cov[!(CNVs_Buko_halleri_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Buko_halleri_HC_GROM_cov$CN_class=="CN0"|CNVs_Buko_halleri_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Buko_halleri_HC_GROM_cov,"Buko_Kowa_halleri_CNVs_overlap_strict.table",row.names=F,sep="\t")
write.xlsx(CNVs_Buko_halleri_HC_GROM_cov,"CNVs_overlap_strict_BukoKowa_halleri.xlsx",sheetName="Buko_Kowa_halleri",row.names=F,overwrite=T)

#merge in Lumpy data
require(tidyr)
BZ_ha_Lumpy_DUP<-read.xlsx("../Lumpy/Buko_Kowa_halleri.vcf_DUP_Lumpy.xlsx",1)
BZ_ha_Lumpy_DEL<-read.xlsx("../Lumpy/Buko_Kowa_halleri.vcf_DEL_Lumpy.xlsx",1)
test1<-separate(data=BZ_ha_Lumpy_DUP,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
BZ_ha_Lumpy_DUP$End<-as.numeric(gsub("END=","",test1$NA4))
test2<-separate(data=BZ_ha_Lumpy_DEL,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
BZ_ha_Lumpy_DEL$End<-as.numeric(gsub("END=","",test2$NA4))
BZ_ha_Lumpy<-rbind(BZ_ha_Lumpy_DUP,BZ_ha_Lumpy_DEL)
BZ_ha_Lumpy<-BZ_ha_Lumpy[(BZ_ha_Lumpy$End-BZ_ha_Lumpy$Pos)<100000,]
BZ_ha_Lumpy$PE_M_std<-BZ_ha_Lumpy$PE_M/(BZ_ha_Lumpy$End-BZ_ha_Lumpy$Pos)
BZ_ha_Lumpy$SU_M_std<-BZ_ha_Lumpy$SU_M/(BZ_ha_Lumpy$End-BZ_ha_Lumpy$Pos)
BZ_ha_Lumpy$SR_M_std<-BZ_ha_Lumpy$SR_M/(BZ_ha_Lumpy$End-BZ_ha_Lumpy$Pos)
BZ_ha_Lumpy$PE_NM_std<-BZ_ha_Lumpy$PE_NM/(BZ_ha_Lumpy$End-BZ_ha_Lumpy$Pos)
BZ_ha_Lumpy$SU_NM_std<-BZ_ha_Lumpy$SU_NM/(BZ_ha_Lumpy$End-BZ_ha_Lumpy$Pos)
BZ_ha_Lumpy$SR_NM_std<-BZ_ha_Lumpy$SR_NM/(BZ_ha_Lumpy$End-BZ_ha_Lumpy$Pos)

BZ_ha_Lumpy_GRange<-GRanges(seqnames=tolower(BZ_ha_Lumpy$Scaffold),ranges=IRanges(start=BZ_ha_Lumpy$Pos,end=BZ_ha_Lumpy$End))
values(BZ_ha_Lumpy_GRange)<-cbind(BZ_ha_Lumpy[,8:ncol(BZ_ha_Lumpy)])

CNVs_Buko_halleri_HC_GROM_cov_GRange<-GRanges(seqnames=tolower(CNVs_Buko_halleri_HC_GROM_cov$Scaffold),ranges=IRanges(start=CNVs_Buko_halleri_HC_GROM_cov$Overlap_start,end=CNVs_Buko_halleri_HC_GROM_cov$Overlap_end))
values(CNVs_Buko_halleri_HC_GROM_cov_GRange)<-cbind(CNVs_Buko_halleri_HC_GROM_cov[,1:54])

Buko_halleri_overlap_strict2<-mergeByOverlaps(CNVs_Buko_halleri_HC_GROM_cov_GRange,BZ_ha_Lumpy_GRange)
Buko_halleri_overlap_strict_df2<-as.data.frame(Buko_halleri_overlap_strict2)[,c(60:113,(113+ncol(BZ_ha_Lumpy)-1):((113+ncol(BZ_ha_Lumpy)-1)+(ncol(BZ_ha_Lumpy)-8)))]

write.xlsx(Buko_halleri_overlap_strict_df2,"CNVs_overlap_strict_BukoKowa_halleri_Lumpy.xlsx",sheetName="Buko_Kowa_halleri",row.names=F)

############################################################
#GROM data#
############################################################
Mias_arenosa_CNVs1<-read.table("Mias.arenosa.cnvs.vcf",sep="\t",header=F)

names(Mias_arenosa_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Mias_arenosa_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Mias_arenosa_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Mias_arenosa_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Mias_arenosa_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Mias_arenosa_CNVs<-data.frame(Mias_arenosa_CNVs1[,1:2],as.numeric(as.character(Mias_arenosa_CNVs_END[,2])),Mias_arenosa_CNVs1[,3:7],Mias_arenosa_CNVs_info)
names(Mias_arenosa_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Mias_arenosa_cnvs_HC_strict<-read.delim("../Cnmops/Kowa/MKar_cnvs_only_strict_minL6_WL350.table",header=T)
Kowa_arenosa_cnvs_HC_strict<-read.delim("../Cnmops/Kowa/KMar_cnvs_only_strict_minL6_WL350.table",header=T)

Mias_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Mias_arenosa_CNVs$Scaffold),ranges=IRanges(start=Mias_arenosa_CNVs$Start_pos,end=Mias_arenosa_CNVs$End_pos))
values(Mias_arenosa_CNVs_GRange)<-Mias_arenosa_CNVs[,4:12]
Kowa_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_arenosa_CNVs$Scaffold),ranges=IRanges(start=Kowa_arenosa_CNVs$Start_pos,end=Kowa_arenosa_CNVs$End_pos))
values(Kowa_arenosa_CNVs_GRange)<-Kowa_arenosa_CNVs[,4:12]
MKar_merge<-mergeByOverlaps(Mias_arenosa_CNVs_GRange,Kowa_arenosa_CNVs_GRange)
MKar_merge_same<-MKar_merge[abs(as.numeric(as.character(MKar_merge@ listData$ Mias_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(MKar_merge@ listData$ Kowa_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
MKar_merge_same_df<-as.data.frame(MKar_merge_same)[,1:3]
names(MKar_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
MKar_merge_same_df_GRange<-GRanges(seqnames=tolower(MKar_merge_same_df$Scaffold),ranges=IRanges(start=MKar_merge_same_df$Start_pos,end=MKar_merge_same_df$End_pos))

#subset Mias_arenosa CNVs to regions without similar CN in Kowa_arenosa

MKar_diff<-setdiff(Mias_arenosa_CNVs_GRange,MKar_merge_same_df_GRange)
MKar<-subsetByOverlaps(Mias_arenosa_CNVs_GRange,MKar_diff)

#overlap cn.mops and GROM regions

Mias_arenosa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Mias_arenosa_cnvs_HC_strict$Chr),ranges=IRanges(start=Mias_arenosa_cnvs_HC_strict$gene_start,end=Mias_arenosa_cnvs_HC_strict$gene_end))
values(Mias_arenosa_cnvs_HC_strict_GRange)<-cbind(Mias_arenosa_cnvs_HC_strict[,1:12],Mias_arenosa_cnvs_HC_strict[,18:40])
Mias_arenosa_overlap_strict<-mergeByOverlaps(Mias_arenosa_cnvs_HC_strict_GRange,MKar)
Mias_arenosa_overlap_strict_df<-as.data.frame(Mias_arenosa_overlap_strict)

Mias_arenosa_CNVs_overlapping_strict<-data.frame(Mias_arenosa_overlap_strict_df[,1:40],Mias_arenosa_overlap_strict_df[,76:89])
names(Mias_arenosa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Mias_arenosa_CNVs_HC_GROM_strict<-rbind(Mias_arenosa_CNVs_overlapping_strict,Kowa_arenosa_CNVs_overlapping_strict)

CNVs_Mias_arenosa_HC_GROM_cov<-Mias_arenosa_CNVs_HC_GROM_strict
CNVs_Mias_arenosa_HC_GROM_cov<-CNVs_Mias_arenosa_HC_GROM_cov[!((CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN3")&CNVs_Mias_arenosa_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Mias_arenosa_HC_GROM_cov<-CNVs_Mias_arenosa_HC_GROM_cov[!(CNVs_Mias_arenosa_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Mias_arenosa_HC_GROM_cov$CN_class=="CN3")),]

write.table(CNVs_Mias_arenosa_HC_GROM_cov,"Mias_Kowa_arenosa_CNVs_overlap_strict.table",row.names=F,sep="\t")
write.xlsx(CNVs_Mias_arenosa_HC_GROM_cov,"CNVs_overlap_strict_MiasKowa_arenosa.xlsx",sheetName="Mias_Kowa_arenosa",row.names=F,overwrite=T)

#merge in Lumpy data
require(tidyr)
MK_ar_Lumpy_DUP<-read.xlsx("../Lumpy/Mias_Kowa_arenosa.vcf_DUP_Lumpy.xlsx",1)
MK_ar_Lumpy_DEL<-read.xlsx("../Lumpy/Mias_Kowa_arenosa.vcf_DEL_Lumpy.xlsx",1)
test1<-separate(data=MK_ar_Lumpy_DUP,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
MK_ar_Lumpy_DUP$End<-as.numeric(gsub("END=","",test1$NA4))
test2<-separate(data=MK_ar_Lumpy_DEL,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
MK_ar_Lumpy_DEL$End<-as.numeric(gsub("END=","",test2$NA4))
MK_ar_Lumpy<-rbind(MK_ar_Lumpy_DUP,MK_ar_Lumpy_DEL)
MK_ar_Lumpy<-MK_ar_Lumpy[(MK_ar_Lumpy$End-MK_ar_Lumpy$Pos)<100000,]
MK_ar_Lumpy$PE_M_std<-MK_ar_Lumpy$PE_M/(MK_ar_Lumpy$End-MK_ar_Lumpy$Pos)
MK_ar_Lumpy$SU_M_std<-MK_ar_Lumpy$SU_M/(MK_ar_Lumpy$End-MK_ar_Lumpy$Pos)
MK_ar_Lumpy$SR_M_std<-MK_ar_Lumpy$SR_M/(MK_ar_Lumpy$End-MK_ar_Lumpy$Pos)
MK_ar_Lumpy$PE_NM_std<-MK_ar_Lumpy$PE_NM/(MK_ar_Lumpy$End-MK_ar_Lumpy$Pos)
MK_ar_Lumpy$SU_NM_std<-MK_ar_Lumpy$SU_NM/(MK_ar_Lumpy$End-MK_ar_Lumpy$Pos)
MK_ar_Lumpy$SR_NM_std<-MK_ar_Lumpy$SR_NM/(MK_ar_Lumpy$End-MK_ar_Lumpy$Pos)

MK_ar_Lumpy_GRange<-GRanges(seqnames=tolower(MK_ar_Lumpy$Scaffold),ranges=IRanges(start=MK_ar_Lumpy$Pos,end=MK_ar_Lumpy$End))
values(MK_ar_Lumpy_GRange)<-cbind(MK_ar_Lumpy[,8:ncol(MK_ar_Lumpy)])

CNVs_Mias_arenosa_HC_GROM_cov_GRange<-GRanges(seqnames=tolower(CNVs_Mias_arenosa_HC_GROM_cov$Scaffold),ranges=IRanges(start=ifelse(CNVs_Mias_arenosa_HC_GROM_cov$Copy_start<CNVs_Mias_arenosa_HC_GROM_cov$GROM_Start_pos,CNVs_Mias_arenosa_HC_GROM_cov$Copy_start,CNVs_Mias_arenosa_HC_GROM_cov$GROM_Start_pos),end=ifelse(CNVs_Mias_arenosa_HC_GROM_cov$Copy_end>CNVs_Mias_arenosa_HC_GROM_cov$GROM_End_pos,CNVs_Mias_arenosa_HC_GROM_cov$Copy_end,CNVs_Mias_arenosa_HC_GROM_cov$GROM_End_pos)))
values(CNVs_Mias_arenosa_HC_GROM_cov_GRange)<-cbind(CNVs_Mias_arenosa_HC_GROM_cov[,1:54])

Mias_arenosa_overlap_strict2<-mergeByOverlaps(CNVs_Mias_arenosa_HC_GROM_cov_GRange,MK_ar_Lumpy_GRange)
Mias_arenosa_overlap_strict_df2<-as.data.frame(Mias_arenosa_overlap_strict2)[,c(60:113,(113+ncol(MK_ar_Lumpy)-1):((113+ncol(MK_ar_Lumpy)-1)+(ncol(MK_ar_Lumpy)-8)))]

write.xlsx(Mias_arenosa_overlap_strict_df2,"CNVs_overlap_strict_MiasKowa_arenosa_Lumpy.xlsx",sheetName="Mias_Kowa_arenosa",row.names=F,overwrite=T)

############################################################
#GROM data#
############################################################
Mias_halleri_CNVs1<-read.table("Mias.halleri.cnvs.vcf",sep="\t",header=F)

names(Mias_halleri_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Mias_halleri_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Mias_halleri_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Mias_halleri_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Mias_halleri_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Mias_halleri_CNVs<-data.frame(Mias_halleri_CNVs1[,1:2],as.numeric(as.character(Mias_halleri_CNVs_END[,2])),Mias_halleri_CNVs1[,3:7],Mias_halleri_CNVs_info)
names(Mias_halleri_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Kowa_halleri_CNVs1<-read.table("Kowa.halleri.cnvs.vcf",sep="\t",header=F)

names(Kowa_halleri_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Kowa_halleri_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_halleri_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Kowa_halleri_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_halleri_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Kowa_halleri_CNVs<-data.frame(Kowa_halleri_CNVs1[,1:2],as.numeric(as.character(Kowa_halleri_CNVs_END[,2])),Kowa_halleri_CNVs1[,3:7],Kowa_halleri_CNVs_info)
names(Kowa_halleri_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Mias_halleri_cnvs_HC_strict<-read.table("../Cnmops/Kowa/MZha_cnvs_only_strict_minL6_WL350.table",header=T)
Kowa_halleri_cnvs_HC_strict<-read.table("../Cnmops/Kowa/ZMha_cnvs_only_strict_minL6_WL350.table",header=T)

Mias_halleri_CNVs_GRange<-GRanges(seqnames=tolower(Mias_halleri_CNVs$Scaffold),ranges=IRanges(start=Mias_halleri_CNVs$Start_pos,end=Mias_halleri_CNVs$End_pos))
values(Mias_halleri_CNVs_GRange)<-Mias_halleri_CNVs[,4:12]
Kowa_halleri_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_halleri_CNVs$Scaffold),ranges=IRanges(start=Kowa_halleri_CNVs$Start_pos,end=Kowa_halleri_CNVs$End_pos))
values(Kowa_halleri_CNVs_GRange)<-Kowa_halleri_CNVs[,4:12]
MZha_merge<-mergeByOverlaps(Mias_halleri_CNVs_GRange,Kowa_halleri_CNVs_GRange)
MZha_merge_same<-MZha_merge[abs(as.numeric(as.character(MZha_merge@ listData$ Mias_halleri_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(MZha_merge@ listData$ Kowa_halleri_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
MZha_merge_same_df<-as.data.frame(MZha_merge_same)[,1:3]
names(MZha_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
MZha_merge_same_df_GRange<-GRanges(seqnames=tolower(MZha_merge_same_df$Scaffold),ranges=IRanges(start=MZha_merge_same_df$Start_pos,end=MZha_merge_same_df$End_pos))

#subset Mias_halleri CNVs to regions without similar CN in Kowa_halleri

MZha_diff<-setdiff(Mias_halleri_CNVs_GRange,MZha_merge_same_df_GRange)
MZha<-subsetByOverlaps(Mias_halleri_CNVs_GRange,MZha_diff)
ZMha_diff<-setdiff(Kowa_halleri_CNVs_GRange,MZha_merge_same_df_GRange)
ZMha<-subsetByOverlaps(Kowa_halleri_CNVs_GRange,ZMha_diff)

#overlap cn.mops and GROM regions

Mias_halleri_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Mias_halleri_cnvs_HC_strict$Chr),ranges=IRanges(start=Mias_halleri_cnvs_HC_strict$gene_start,end=Mias_halleri_cnvs_HC_strict$gene_end))
values(Mias_halleri_cnvs_HC_strict_GRange)<-cbind(Mias_halleri_cnvs_HC_strict[,1:12],Mias_halleri_cnvs_HC_strict[,18:40])
Mias_halleri_overlap_strict<-mergeByOverlaps(Mias_halleri_cnvs_HC_strict_GRange,MZha)
Mias_halleri_overlap_strict_df<-as.data.frame(Mias_halleri_overlap_strict)
Kowa_halleri_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Kowa_halleri_cnvs_HC_strict$Chr),ranges=IRanges(start=Kowa_halleri_cnvs_HC_strict$gene_start,end=Kowa_halleri_cnvs_HC_strict$gene_end))
values(Kowa_halleri_cnvs_HC_strict_GRange)<-cbind(Kowa_halleri_cnvs_HC_strict[,1:12],Kowa_halleri_cnvs_HC_strict[,18:40])
Kowa_halleri_overlap_strict<-mergeByOverlaps(Kowa_halleri_cnvs_HC_strict_GRange,ZMha)
Kowa_halleri_overlap_strict_df<-as.data.frame(Kowa_halleri_overlap_strict)

Mias_halleri_CNVs_overlapping_strict<-data.frame(Mias_halleri_overlap_strict_df[,1:40],Mias_halleri_overlap_strict_df[,76:89])
names(Mias_halleri_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")
Kowa_halleri_CNVs_overlapping_strict<-data.frame(Kowa_halleri_overlap_strict_df[,1:40],Kowa_halleri_overlap_strict_df[,76:89])
names(Kowa_halleri_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Mias_halleri_CNVs_HC_GROM_strict<-rbind(Mias_halleri_CNVs_overlapping_strict,Kowa_halleri_CNVs_overlapping_strict)

CNVs_Mias_halleri_HC_GROM_cov<-Mias_halleri_CNVs_HC_GROM_strict
CNVs_Mias_halleri_HC_GROM_cov<-CNVs_Mias_halleri_HC_GROM_cov[!((CNVs_Mias_halleri_HC_GROM_cov$CN_class=="CN0"|CNVs_Mias_halleri_HC_GROM_cov$CN_class=="CN1")&CNVs_Mias_halleri_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Mias_halleri_HC_GROM_cov<-CNVs_Mias_halleri_HC_GROM_cov[!(CNVs_Mias_halleri_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Mias_halleri_HC_GROM_cov$CN_class=="CN0"|CNVs_Mias_halleri_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Mias_halleri_HC_GROM_cov,"Mias_Kowa_halleri_CNVs_overlap_strict.table",row.names=F,sep="\t")
write.xlsx(CNVs_Mias_halleri_HC_GROM_cov,"CNVs_overlap_strict_MiasKowa_halleri.xlsx",sheetName="Mias_Kowa_halleri",row.names=F,overwrite=T)

#merge in Lumpy data
require(tidyr)
MZ_ha_Lumpy_DUP<-read.xlsx("../Lumpy/Mias_Kowa_halleri.vcf_DUP_Lumpy.xlsx",1)
MZ_ha_Lumpy_DEL<-read.xlsx("../Lumpy/Mias_Kowa_halleri.vcf_DEL_Lumpy.xlsx",1)
test1<-separate(data=MZ_ha_Lumpy_DUP,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
MZ_ha_Lumpy_DUP$End<-as.numeric(gsub("END=","",test1$NA4))
test2<-separate(data=MZ_ha_Lumpy_DEL,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
MZ_ha_Lumpy_DEL$End<-as.numeric(gsub("END=","",test2$NA4))
MZ_ha_Lumpy<-rbind(MZ_ha_Lumpy_DUP,MZ_ha_Lumpy_DEL)
MZ_ha_Lumpy<-MZ_ha_Lumpy[(MZ_ha_Lumpy$End-MZ_ha_Lumpy$Pos)<100000,]
MZ_ha_Lumpy$PE_M_std<-MZ_ha_Lumpy$PE_M/(MZ_ha_Lumpy$End-MZ_ha_Lumpy$Pos)
MZ_ha_Lumpy$SU_M_std<-MZ_ha_Lumpy$SU_M/(MZ_ha_Lumpy$End-MZ_ha_Lumpy$Pos)
MZ_ha_Lumpy$SR_M_std<-MZ_ha_Lumpy$SR_M/(MZ_ha_Lumpy$End-MZ_ha_Lumpy$Pos)
MZ_ha_Lumpy$PE_NM_std<-MZ_ha_Lumpy$PE_NM/(MZ_ha_Lumpy$End-MZ_ha_Lumpy$Pos)
MZ_ha_Lumpy$SU_NM_std<-MZ_ha_Lumpy$SU_NM/(MZ_ha_Lumpy$End-MZ_ha_Lumpy$Pos)
MZ_ha_Lumpy$SR_NM_std<-MZ_ha_Lumpy$SR_NM/(MZ_ha_Lumpy$End-MZ_ha_Lumpy$Pos)

MZ_ha_Lumpy_GRange<-GRanges(seqnames=tolower(MZ_ha_Lumpy$Scaffold),ranges=IRanges(start=MZ_ha_Lumpy$Pos,end=MZ_ha_Lumpy$End))
values(MZ_ha_Lumpy_GRange)<-cbind(MZ_ha_Lumpy[,8:ncol(MZ_ha_Lumpy)])

CNVs_Mias_halleri_HC_GROM_cov_GRange<-GRanges(seqnames=tolower(CNVs_Mias_halleri_HC_GROM_cov$Scaffold),ranges=IRanges(start=CNVs_Mias_halleri_HC_GROM_cov$Overlap_start,end=CNVs_Mias_halleri_HC_GROM_cov$Overlap_end))
values(CNVs_Mias_halleri_HC_GROM_cov_GRange)<-cbind(CNVs_Mias_halleri_HC_GROM_cov[,1:54])

Mias_halleri_overlap_strict2<-mergeByOverlaps(CNVs_Mias_halleri_HC_GROM_cov_GRange,MZ_ha_Lumpy_GRange)
Mias_halleri_overlap_strict_df2<-as.data.frame(Mias_halleri_overlap_strict2)[,c(60:113,(113+ncol(MZ_ha_Lumpy)-1):((113+ncol(MZ_ha_Lumpy)-1)+(ncol(MZ_ha_Lumpy)-8)))]

write.xlsx(Mias_halleri_overlap_strict_df2,"CNVs_overlap_strict_MiasKowa_halleri_Lumpy.xlsx",sheetName="Mias_Kowa_halleri",row.names=F)

############################################################
#GROM data#
############################################################
Piek_arenosa_CNVs1<-read.table("Piek.arenosa.cnvs.vcf",sep="\t",header=F)

names(Piek_arenosa_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Piek_arenosa_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Piek_arenosa_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Piek_arenosa_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Piek_arenosa_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Piek_arenosa_CNVs<-data.frame(Piek_arenosa_CNVs1[,1:2],as.numeric(as.character(Piek_arenosa_CNVs_END[,2])),Piek_arenosa_CNVs1[,3:7],Piek_arenosa_CNVs_info)
names(Piek_arenosa_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Piek_arenosa_cnvs_HC_strict<-read.table("../Cnmops/Kowa/PKar_cnvs_only_strict_minL6_WL350.table",header=T)
Kowa_arenosa_cnvs_HC_strict<-read.table("../Cnmops/Kowa/KPar_cnvs_only_strict_minL6_WL350.table",header=T)

require(GenomicRanges)

Piek_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Piek_arenosa_CNVs$Scaffold),ranges=IRanges(start=Piek_arenosa_CNVs$Start_pos,end=Piek_arenosa_CNVs$End_pos))
values(Piek_arenosa_CNVs_GRange)<-Piek_arenosa_CNVs[,4:12]
Kowa_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_arenosa_CNVs$Scaffold),ranges=IRanges(start=Kowa_arenosa_CNVs$Start_pos,end=Kowa_arenosa_CNVs$End_pos))
values(Kowa_arenosa_CNVs_GRange)<-Kowa_arenosa_CNVs[,4:12]
PKar_merge<-mergeByOverlaps(Piek_arenosa_CNVs_GRange,Kowa_arenosa_CNVs_GRange)
PKar_merge_same<-PKar_merge[abs(as.numeric(as.character(PKar_merge@ listData$ Piek_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(PKar_merge@ listData$ Kowa_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
PKar_merge_same_df<-as.data.frame(PKar_merge_same)[,1:3]
names(PKar_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
PKar_merge_same_df_GRange<-GRanges(seqnames=tolower(PKar_merge_same_df$Scaffold),ranges=IRanges(start=PKar_merge_same_df$Start_pos,end=PKar_merge_same_df$End_pos))

#subset Piek_arenosa CNVs to regions without similar CN in Kowa_arenosa

PKar_diff<-setdiff(Piek_arenosa_CNVs_GRange,PKar_merge_same_df_GRange)
PKar<-subsetByOverlaps(Piek_arenosa_CNVs_GRange,PKar_diff)
KPar_diff<-setdiff(Kowa_arenosa_CNVs_GRange,PKar_merge_same_df_GRange)
KPar<-subsetByOverlaps(Kowa_arenosa_CNVs_GRange,KPar_diff)

#overlap cn.mops and GROM regions

Piek_arenosa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Piek_arenosa_cnvs_HC_strict$Chr),ranges=IRanges(start=Piek_arenosa_cnvs_HC_strict$gene_start,end=Piek_arenosa_cnvs_HC_strict$gene_end))
values(Piek_arenosa_cnvs_HC_strict_GRange)<-cbind(Piek_arenosa_cnvs_HC_strict[,1:12],Piek_arenosa_cnvs_HC_strict[,18:40])
Piek_arenosa_overlap_strict<-mergeByOverlaps(Piek_arenosa_cnvs_HC_strict_GRange,PKar)
Piek_arenosa_overlap_strict_df<-as.data.frame(Piek_arenosa_overlap_strict)
Kowa_arenosa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Kowa_arenosa_cnvs_HC_strict$Chr),ranges=IRanges(start=Kowa_arenosa_cnvs_HC_strict$gene_start,end=Kowa_arenosa_cnvs_HC_strict$gene_end))
values(Kowa_arenosa_cnvs_HC_strict_GRange)<-cbind(Kowa_arenosa_cnvs_HC_strict[,1:12],Kowa_arenosa_cnvs_HC_strict[,18:40])
Kowa_arenosa_overlap_strict<-mergeByOverlaps(Kowa_arenosa_cnvs_HC_strict_GRange,KPar)
Kowa_arenosa_overlap_strict_df<-as.data.frame(Kowa_arenosa_overlap_strict)

Piek_arenosa_CNVs_overlapping_strict<-data.frame(Piek_arenosa_overlap_strict_df[,1:40],Piek_arenosa_overlap_strict_df[,76:89])
names(Piek_arenosa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")
Kowa_arenosa_CNVs_overlapping_strict<-data.frame(Kowa_arenosa_overlap_strict_df[,1:40],Kowa_arenosa_overlap_strict_df[,76:89])
names(Kowa_arenosa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Piek_arenosa_CNVs_HC_GROM_strict<-rbind(Piek_arenosa_CNVs_overlapping_strict,Kowa_arenosa_CNVs_overlapping_strict)

CNVs_Piek_arenosa_HC_GROM_cov<-Piek_arenosa_CNVs_HC_GROM_strict
CNVs_Piek_arenosa_HC_GROM_cov<-CNVs_Piek_arenosa_HC_GROM_cov[!((CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN3")&CNVs_Piek_arenosa_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Piek_arenosa_HC_GROM_cov<-CNVs_Piek_arenosa_HC_GROM_cov[!(CNVs_Piek_arenosa_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Piek_arenosa_HC_GROM_cov$CN_class=="CN3")),]

write.table(CNVs_Piek_arenosa_HC_GROM_cov,"Piek_Kowa_arenosa_CNVs_overlap_strict.table",row.names=F,sep="\t")
write.xlsx(CNVs_Piek_arenosa_HC_GROM_cov,"CNVs_overlap_strict_PiekKowa_arenosa.xlsx",sheetName="Piek_Kowa_arenosa",row.names=F,overwrite=T)

#merge in Lumpy data
require(tidyr)
PK_ar_Lumpy_DUP<-read.xlsx("../Lumpy/Piek_Kowa_arenosa.vcf_DUP_Lumpy.xlsx",1)
PK_ar_Lumpy_DEL<-read.xlsx("../Lumpy/Piek_Kowa_arenosa.vcf_DEL_Lumpy.xlsx",1)
test1<-separate(data=PK_ar_Lumpy_DUP,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
PK_ar_Lumpy_DUP$End<-as.numeric(gsub("END=","",test1$NA4))
test2<-separate(data=PK_ar_Lumpy_DEL,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
PK_ar_Lumpy_DEL$End<-as.numeric(gsub("END=","",test2$NA4))
PK_ar_Lumpy<-rbind(PK_ar_Lumpy_DUP,PK_ar_Lumpy_DEL)
PK_ar_Lumpy<-PK_ar_Lumpy[(PK_ar_Lumpy$End-PK_ar_Lumpy$Pos)<100000,]
PK_ar_Lumpy$PE_M_std<-PK_ar_Lumpy$PE_M/(PK_ar_Lumpy$End-PK_ar_Lumpy$Pos)
PK_ar_Lumpy$SU_M_std<-PK_ar_Lumpy$SU_M/(PK_ar_Lumpy$End-PK_ar_Lumpy$Pos)
PK_ar_Lumpy$SR_M_std<-PK_ar_Lumpy$SR_M/(PK_ar_Lumpy$End-PK_ar_Lumpy$Pos)
PK_ar_Lumpy$PE_NM_std<-PK_ar_Lumpy$PE_NM/(PK_ar_Lumpy$End-PK_ar_Lumpy$Pos)
PK_ar_Lumpy$SU_NM_std<-PK_ar_Lumpy$SU_NM/(PK_ar_Lumpy$End-PK_ar_Lumpy$Pos)
PK_ar_Lumpy$SR_NM_std<-PK_ar_Lumpy$SR_NM/(PK_ar_Lumpy$End-PK_ar_Lumpy$Pos)

PK_ar_Lumpy_GRange<-GRanges(seqnames=tolower(PK_ar_Lumpy$Scaffold),ranges=IRanges(start=PK_ar_Lumpy$Pos,end=PK_ar_Lumpy$End))
values(PK_ar_Lumpy_GRange)<-cbind(PK_ar_Lumpy[,8:ncol(PK_ar_Lumpy)])

CNVs_Piek_arenosa_HC_GROM_cov_GRange<-GRanges(seqnames=tolower(CNVs_Piek_arenosa_HC_GROM_cov$Scaffold),ranges=IRanges(start=ifelse(CNVs_Piek_arenosa_HC_GROM_cov$Copy_start<CNVs_Piek_arenosa_HC_GROM_cov$GROM_Start_pos,CNVs_Piek_arenosa_HC_GROM_cov$Copy_start,CNVs_Piek_arenosa_HC_GROM_cov$GROM_Start_pos),end=ifelse(CNVs_Piek_arenosa_HC_GROM_cov$Copy_end>CNVs_Piek_arenosa_HC_GROM_cov$GROM_End_pos,CNVs_Piek_arenosa_HC_GROM_cov$Copy_end,CNVs_Piek_arenosa_HC_GROM_cov$GROM_End_pos)))
values(CNVs_Piek_arenosa_HC_GROM_cov_GRange)<-cbind(CNVs_Piek_arenosa_HC_GROM_cov[,1:54])

Piek_arenosa_overlap_strict2<-mergeByOverlaps(CNVs_Piek_arenosa_HC_GROM_cov_GRange,PK_ar_Lumpy_GRange)
Piek_arenosa_overlap_strict_df2<-as.data.frame(Piek_arenosa_overlap_strict2)[,c(60:113,(113+ncol(PK_ar_Lumpy)-1):((113+ncol(PK_ar_Lumpy)-1)+(ncol(PK_ar_Lumpy)-8)))]

write.xlsx(Piek_arenosa_overlap_strict_df2,"CNVs_overlap_strict_PiekKowa_arenosa_Lumpy.xlsx",sheetName="Piek_Kowa_arenosa",row.names=F,overwrite=T)

############################################################
#GROM data#
############################################################
Piek_halleri_CNVs1<-read.table("Piek.halleri.cnvs.vcf",sep="\t",header=F)

names(Piek_halleri_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Piek_halleri_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Piek_halleri_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Piek_halleri_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Piek_halleri_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Piek_halleri_CNVs<-data.frame(Piek_halleri_CNVs1[,1:2],as.numeric(as.character(Piek_halleri_CNVs_END[,2])),Piek_halleri_CNVs1[,3:7],Piek_halleri_CNVs_info)
names(Piek_halleri_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Kowa_halleri_CNVs1<-read.table("Kowa.halleri.cnvs.vcf",sep="\t",header=F)

names(Kowa_halleri_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Kowa_halleri_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_halleri_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Kowa_halleri_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_halleri_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Kowa_halleri_CNVs<-data.frame(Kowa_halleri_CNVs1[,1:2],as.numeric(as.character(Kowa_halleri_CNVs_END[,2])),Kowa_halleri_CNVs1[,3:7],Kowa_halleri_CNVs_info)
names(Kowa_halleri_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Piek_halleri_cnvs_HC_strict<-read.table("../Cnmops/Kowa/PZha_cnvs_only_strict_minL6_WL350.table",header=T)
Kowa_halleri_cnvs_HC_strict<-read.table("../Cnmops/Kowa/ZPha_cnvs_only_strict_minL6_WL350.table",header=T)

Piek_halleri_CNVs_GRange<-GRanges(seqnames=tolower(Piek_halleri_CNVs$Scaffold),ranges=IRanges(start=Piek_halleri_CNVs$Start_pos,end=Piek_halleri_CNVs$End_pos))
values(Piek_halleri_CNVs_GRange)<-Piek_halleri_CNVs[,4:12]
Kowa_halleri_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_halleri_CNVs$Scaffold),ranges=IRanges(start=Kowa_halleri_CNVs$Start_pos,end=Kowa_halleri_CNVs$End_pos))
values(Kowa_halleri_CNVs_GRange)<-Kowa_halleri_CNVs[,4:12]
PZha_merge<-mergeByOverlaps(Piek_halleri_CNVs_GRange,Kowa_halleri_CNVs_GRange)
PZha_merge_same<-PZha_merge[abs(as.numeric(as.character(PZha_merge@ listData$ Piek_halleri_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(PZha_merge@ listData$ Kowa_halleri_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
PZha_merge_same_df<-as.data.frame(PZha_merge_same)[,1:3]
names(PZha_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
PZha_merge_same_df_GRange<-GRanges(seqnames=tolower(PZha_merge_same_df$Scaffold),ranges=IRanges(start=PZha_merge_same_df$Start_pos,end=PZha_merge_same_df$End_pos))

#subset Piek_halleri CNVs to regions without similar CN in Kowa_halleri

PZha_diff<-setdiff(Piek_halleri_CNVs_GRange,PZha_merge_same_df_GRange)
PZha<-subsetByOverlaps(Piek_halleri_CNVs_GRange,PZha_diff)
ZPha_diff<-setdiff(Kowa_halleri_CNVs_GRange,PZha_merge_same_df_GRange)
ZPha<-subsetByOverlaps(Kowa_halleri_CNVs_GRange,ZPha_diff)

#overlap cn.mops and GROM regions

Piek_halleri_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Piek_halleri_cnvs_HC_strict$Chr),ranges=IRanges(start=Piek_halleri_cnvs_HC_strict$gene_start,end=Piek_halleri_cnvs_HC_strict$gene_end))
values(Piek_halleri_cnvs_HC_strict_GRange)<-cbind(Piek_halleri_cnvs_HC_strict[,1:12],Piek_halleri_cnvs_HC_strict[,18:40])
Piek_halleri_overlap_strict<-mergeByOverlaps(Piek_halleri_cnvs_HC_strict_GRange,PZha)
Piek_halleri_overlap_strict_df<-as.data.frame(Piek_halleri_overlap_strict)
Kowa_halleri_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Kowa_halleri_cnvs_HC_strict$Chr),ranges=IRanges(start=Kowa_halleri_cnvs_HC_strict$gene_start,end=Kowa_halleri_cnvs_HC_strict$gene_end))
values(Kowa_halleri_cnvs_HC_strict_GRange)<-cbind(Kowa_halleri_cnvs_HC_strict[,1:12],Kowa_halleri_cnvs_HC_strict[,18:40])
Kowa_halleri_overlap_strict<-mergeByOverlaps(Kowa_halleri_cnvs_HC_strict_GRange,ZPha)
Kowa_halleri_overlap_strict_df<-as.data.frame(Kowa_halleri_overlap_strict)

Piek_halleri_CNVs_overlapping_strict<-data.frame(Piek_halleri_overlap_strict_df[,1:40],Piek_halleri_overlap_strict_df[,76:89])
names(Piek_halleri_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")
Kowa_halleri_CNVs_overlapping_strict<-data.frame(Kowa_halleri_overlap_strict_df[,1:40],Kowa_halleri_overlap_strict_df[,76:89])
names(Kowa_halleri_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Piek_halleri_CNVs_HC_GROM_strict<-rbind(Piek_halleri_CNVs_overlapping_strict,Kowa_halleri_CNVs_overlapping_strict)

CNVs_Piek_halleri_HC_GROM_cov<-Piek_halleri_CNVs_HC_GROM_strict
CNVs_Piek_halleri_HC_GROM_cov<-CNVs_Piek_halleri_HC_GROM_cov[!((CNVs_Piek_halleri_HC_GROM_cov$CN_class=="CN0"|CNVs_Piek_halleri_HC_GROM_cov$CN_class=="CN1")&CNVs_Piek_halleri_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Piek_halleri_HC_GROM_cov<-CNVs_Piek_halleri_HC_GROM_cov[!(CNVs_Piek_halleri_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Piek_halleri_HC_GROM_cov$CN_class=="CN0"|CNVs_Piek_halleri_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Piek_halleri_HC_GROM_cov,"Piek_Kowa_halleri_CNVs_overlap_strict.table",row.names=F,sep="\t")
write.xlsx(CNVs_Piek_halleri_HC_GROM_cov,"CNVs_overlap_strict_PiekKowa_halleri.xlsx",sheetName="Piek_Kowa_halleri",row.names=F,overwrite=T)

#merge in Lumpy data
require(tidyr)
PZ_ha_Lumpy_DUP<-read.xlsx("../Lumpy/Piek_Kowa_halleri.vcf_DUP_Lumpy.xlsx",1)
PZ_ha_Lumpy_DEL<-read.xlsx("../Lumpy/Piek_Kowa_halleri.vcf_DEL_Lumpy.xlsx",1)
test1<-separate(data=PZ_ha_Lumpy_DUP,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
PZ_ha_Lumpy_DUP$End<-as.numeric(gsub("END=","",test1$NA4))
test2<-separate(data=PZ_ha_Lumpy_DEL,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
PZ_ha_Lumpy_DEL$End<-as.numeric(gsub("END=","",test2$NA4))
PZ_ha_Lumpy<-rbind(PZ_ha_Lumpy_DUP,PZ_ha_Lumpy_DEL)
PZ_ha_Lumpy<-PZ_ha_Lumpy[(PZ_ha_Lumpy$End-PZ_ha_Lumpy$Pos)<100000,]
PZ_ha_Lumpy$PE_M_std<-PZ_ha_Lumpy$PE_M/(PZ_ha_Lumpy$End-PZ_ha_Lumpy$Pos)
PZ_ha_Lumpy$SU_M_std<-PZ_ha_Lumpy$SU_M/(PZ_ha_Lumpy$End-PZ_ha_Lumpy$Pos)
PZ_ha_Lumpy$SR_M_std<-PZ_ha_Lumpy$SR_M/(PZ_ha_Lumpy$End-PZ_ha_Lumpy$Pos)
PZ_ha_Lumpy$PE_NM_std<-PZ_ha_Lumpy$PE_NM/(PZ_ha_Lumpy$End-PZ_ha_Lumpy$Pos)
PZ_ha_Lumpy$SU_NM_std<-PZ_ha_Lumpy$SU_NM/(PZ_ha_Lumpy$End-PZ_ha_Lumpy$Pos)
PZ_ha_Lumpy$SR_NM_std<-PZ_ha_Lumpy$SR_NM/(PZ_ha_Lumpy$End-PZ_ha_Lumpy$Pos)

PZ_ha_Lumpy_GRange<-GRanges(seqnames=tolower(PZ_ha_Lumpy$Scaffold),ranges=IRanges(start=PZ_ha_Lumpy$Pos,end=PZ_ha_Lumpy$End))
values(PZ_ha_Lumpy_GRange)<-cbind(PZ_ha_Lumpy[,8:ncol(PZ_ha_Lumpy)])

CNVs_Piek_halleri_HC_GROM_cov_GRange<-GRanges(seqnames=tolower(CNVs_Piek_halleri_HC_GROM_cov$Scaffold),ranges=IRanges(start=CNVs_Piek_halleri_HC_GROM_cov$Overlap_start,end=CNVs_Piek_halleri_HC_GROM_cov$Overlap_end))
values(CNVs_Piek_halleri_HC_GROM_cov_GRange)<-cbind(CNVs_Piek_halleri_HC_GROM_cov[,1:54])

Piek_halleri_overlap_strict2<-mergeByOverlaps(CNVs_Piek_halleri_HC_GROM_cov_GRange,PZ_ha_Lumpy_GRange)
Piek_halleri_overlap_strict_df2<-as.data.frame(Piek_halleri_overlap_strict2)[,c(60:113,(113+ncol(PZ_ha_Lumpy)-1):((113+ncol(PZ_ha_Lumpy)-1)+(ncol(PZ_ha_Lumpy)-8)))]

write.xlsx(Piek_halleri_overlap_strict_df2,"CNVs_overlap_strict_PiekKowa_halleri_Lumpy.xlsx",sheetName="Piek_Kowa_halleri",row.names=F)

############################################################
#GROM data#
############################################################
Kato_arenosa_CNVs1<-read.table("Kato.arenosa.cnvs.vcf",sep="\t",header=F)

names(Kato_arenosa_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Kato_arenosa_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Kato_arenosa_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Kato_arenosa_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Kato_arenosa_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Kato_arenosa_CNVs<-data.frame(Kato_arenosa_CNVs1[,1:2],as.numeric(as.character(Kato_arenosa_CNVs_END[,2])),Kato_arenosa_CNVs1[,3:7],Kato_arenosa_CNVs_info)
names(Kato_arenosa_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Kato_arenosa_cnvs_HC_strict<-read.table("../Cnmops/Kowa/KaKoar_cnvs_only_strict_minL6_WL350.table",header=T)
Kowa_arenosa_cnvs_HC_strict<-read.table("../Cnmops/Kowa/KoKaar_cnvs_only_strict_minL6_WL350.table",header=T)

Kato_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Kato_arenosa_CNVs$Scaffold),ranges=IRanges(start=Kato_arenosa_CNVs$Start_pos,end=Kato_arenosa_CNVs$End_pos))
values(Kato_arenosa_CNVs_GRange)<-Kato_arenosa_CNVs[,4:12]
Kowa_arenosa_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_arenosa_CNVs$Scaffold),ranges=IRanges(start=Kowa_arenosa_CNVs$Start_pos,end=Kowa_arenosa_CNVs$End_pos))
values(Kowa_arenosa_CNVs_GRange)<-Kowa_arenosa_CNVs[,4:12]
KaKoar_merge<-mergeByOverlaps(Kato_arenosa_CNVs_GRange,Kowa_arenosa_CNVs_GRange)
KaKoar_merge_same<-KaKoar_merge[abs(as.numeric(as.character(KaKoar_merge@ listData$ Kato_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(KaKoar_merge@ listData$ Kowa_arenosa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
KaKoar_merge_same_df<-as.data.frame(KaKoar_merge_same)[,1:3]
names(KaKoar_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
KaKoar_merge_same_df_GRange<-GRanges(seqnames=tolower(KaKoar_merge_same_df$Scaffold),ranges=IRanges(start=KaKoar_merge_same_df$Start_pos,end=KaKoar_merge_same_df$End_pos))

#subset Kato_arenosa CNVs to regions without similar CN in Kowa_arenosa

KaKoar_diff<-setdiff(Kato_arenosa_CNVs_GRange,KaKoar_merge_same_df_GRange)
KaKoar<-subsetByOverlaps(Kato_arenosa_CNVs_GRange,KaKoar_diff)
KoKaar_diff<-setdiff(Kowa_arenosa_CNVs_GRange,KaKoar_merge_same_df_GRange)
KoKaar<-subsetByOverlaps(Kowa_arenosa_CNVs_GRange,KoKaar_diff)

#overlap cn.mops and GROM regions

Kato_arenosa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Kato_arenosa_cnvs_HC_strict$Chr),ranges=IRanges(start=Kato_arenosa_cnvs_HC_strict$gene_start,end=Kato_arenosa_cnvs_HC_strict$gene_end))
values(Kato_arenosa_cnvs_HC_strict_GRange)<-cbind(Kato_arenosa_cnvs_HC_strict[,1:12],Kato_arenosa_cnvs_HC_strict[,18:40])
Kato_arenosa_overlap_strict<-mergeByOverlaps(Kato_arenosa_cnvs_HC_strict_GRange,KaKoar)
Kato_arenosa_overlap_strict_df<-as.data.frame(Kato_arenosa_overlap_strict)
Kowa_arenosa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Kowa_arenosa_cnvs_HC_strict$Chr),ranges=IRanges(start=Kowa_arenosa_cnvs_HC_strict$gene_start,end=Kowa_arenosa_cnvs_HC_strict$gene_end))
values(Kowa_arenosa_cnvs_HC_strict_GRange)<-cbind(Kowa_arenosa_cnvs_HC_strict[,1:12],Kowa_arenosa_cnvs_HC_strict[,18:40])
Kowa_arenosa_overlap_strict<-mergeByOverlaps(Kowa_arenosa_cnvs_HC_strict_GRange,KoKaar)
Kowa_arenosa_overlap_strict_df<-as.data.frame(Kowa_arenosa_overlap_strict)

Kato_arenosa_CNVs_overlapping_strict<-data.frame(Kato_arenosa_overlap_strict_df[,1:40],Kato_arenosa_overlap_strict_df[,76:89])
names(Kato_arenosa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")
Kowa_arenosa_CNVs_overlapping_strict<-data.frame(Kowa_arenosa_overlap_strict_df[,1:40],Kowa_arenosa_overlap_strict_df[,76:89])
names(Kowa_arenosa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Lyr_Gene","Copy_scaffold","Copy_start","Copy_end","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Widthofoverlap","Ath_ID","Type","Short_description","Curator_summary","Computational_description","Araport11_short_name","GO","Mapman_category","Orthofinder","Name","Function","Localisation","EC.TC","DOI","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Kato_arenosa_CNVs_HC_GROM_strict<-rbind(Kato_arenosa_CNVs_overlapping_strict,Kowa_arenosa_CNVs_overlapping_strict)

CNVs_Kato_arenosa_HC_GROM_cov<-Kato_arenosa_CNVs_HC_GROM_strict
CNVs_Kato_arenosa_HC_GROM_cov<-CNVs_Kato_arenosa_HC_GROM_cov[!((CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN3")&CNVs_Kato_arenosa_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Kato_arenosa_HC_GROM_cov<-CNVs_Kato_arenosa_HC_GROM_cov[!(CNVs_Kato_arenosa_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN0"|CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN1"|CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN2"|CNVs_Kato_arenosa_HC_GROM_cov$CN_class=="CN3")),]

write.table(CNVs_Kato_arenosa_HC_GROM_cov,"Kato_Kowa_arenosa_CNVs_overlap_strict.table",row.names=F,sep="\t")
write.xlsx(CNVs_Kato_arenosa_HC_GROM_cov,"CNVs_overlap_strict_KatoKowa_arenosa.xlsx",sheetName="Kato_Kowa_arenosa",row.names=F,overwrite=T)

#merge in Lumpy data
require(tidyr)
KaKo_ar_Lumpy_DUP<-read.xlsx("../Lumpy/Kato_Kowa_arenosa.vcf_DUP_Lumpy.xlsx",1)
KaKo_ar_Lumpy_DEL<-read.xlsx("../Lumpy/Kato_Kowa_arenosa.vcf_DEL_Lumpy.xlsx",1)
test1<-separate(data=KaKo_ar_Lumpy_DUP,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
KaKo_ar_Lumpy_DUP$End<-as.numeric(gsub("END=","",test1$NA4))
test2<-separate(data=KaKo_ar_Lumpy_DEL,into=c(paste0(rep("NA",12),c(1:12))),col=INFO,sep=";",fill="right")
KaKo_ar_Lumpy_DEL$End<-as.numeric(gsub("END=","",test2$NA4))
KaKo_ar_Lumpy<-rbind(KaKo_ar_Lumpy_DUP,KaKo_ar_Lumpy_DEL)
KaKo_ar_Lumpy<-KaKo_ar_Lumpy[(KaKo_ar_Lumpy$End-KaKo_ar_Lumpy$Pos)<100000,]
KaKo_ar_Lumpy$PE_M_std<-KaKo_ar_Lumpy$PE_M/(KaKo_ar_Lumpy$End-KaKo_ar_Lumpy$Pos)
KaKo_ar_Lumpy$SU_M_std<-KaKo_ar_Lumpy$SU_M/(KaKo_ar_Lumpy$End-KaKo_ar_Lumpy$Pos)
KaKo_ar_Lumpy$SR_M_std<-KaKo_ar_Lumpy$SR_M/(KaKo_ar_Lumpy$End-KaKo_ar_Lumpy$Pos)
KaKo_ar_Lumpy$PE_NM_std<-KaKo_ar_Lumpy$PE_NM/(KaKo_ar_Lumpy$End-KaKo_ar_Lumpy$Pos)
KaKo_ar_Lumpy$SU_NM_std<-KaKo_ar_Lumpy$SU_NM/(KaKo_ar_Lumpy$End-KaKo_ar_Lumpy$Pos)
KaKo_ar_Lumpy$SR_NM_std<-KaKo_ar_Lumpy$SR_NM/(KaKo_ar_Lumpy$End-KaKo_ar_Lumpy$Pos)

KaKo_ar_Lumpy_GRange<-GRanges(seqnames=tolower(KaKo_ar_Lumpy$Scaffold),ranges=IRanges(start=KaKo_ar_Lumpy$Pos,end=KaKo_ar_Lumpy$End))
values(KaKo_ar_Lumpy_GRange)<-cbind(KaKo_ar_Lumpy[,8:ncol(KaKo_ar_Lumpy)])

CNVs_Kato_arenosa_HC_GROM_cov_GRange<-GRanges(seqnames=tolower(CNVs_Kato_arenosa_HC_GROM_cov$Scaffold),ranges=IRanges(start=ifelse(CNVs_Kato_arenosa_HC_GROM_cov$Copy_start<CNVs_Kato_arenosa_HC_GROM_cov$GROM_Start_pos,CNVs_Kato_arenosa_HC_GROM_cov$Copy_start,CNVs_Kato_arenosa_HC_GROM_cov$GROM_Start_pos),end=ifelse(CNVs_Kato_arenosa_HC_GROM_cov$Copy_end>CNVs_Kato_arenosa_HC_GROM_cov$GROM_End_pos,CNVs_Kato_arenosa_HC_GROM_cov$Copy_end,CNVs_Kato_arenosa_HC_GROM_cov$GROM_End_pos)))
values(CNVs_Kato_arenosa_HC_GROM_cov_GRange)<-cbind(CNVs_Kato_arenosa_HC_GROM_cov[,1:54])

Kato_arenosa_overlap_strict2<-mergeByOverlaps(CNVs_Kato_arenosa_HC_GROM_cov_GRange,KaKo_ar_Lumpy_GRange)
Kato_arenosa_overlap_strict_df2<-as.data.frame(Kato_arenosa_overlap_strict2)[,c(60:113,(113+ncol(KaKo_ar_Lumpy)-1):((113+ncol(KaKo_ar_Lumpy)-1)+(ncol(KaKo_ar_Lumpy)-8)))]

write.xlsx(Kato_arenosa_overlap_strict_df2,"CNVs_overlap_strict_KatoKowa_arenosa_Lumpy.xlsx",sheetName="Kato_Kowa_arenosa",row.names=F,overwrite=T)

#exchange for new orthologues

Desc_HM<-merge(Thal_description,UtesMetal,by.x="Ath_ID",by.y="AGI_Number",all.x=T)
Kato_arenosa_overlap_strict_df3<-Kato_arenosa_overlap_strict_df2[,-c(19:40)]
Kato_arenosa_overlap_strict_df4<-merge(Desc_HM,Kato_arenosa_overlap_strict_df3,by.x="Alyr_ID",by.y="Lyr_Gene",all.y=T)

Buko_arenosa_overlap_strict_df3<-Buko_arenosa_overlap_strict_df2[,-c(19:40)]
Buko_arenosa_overlap_strict_df4<-merge(Desc_HM,Buko_arenosa_overlap_strict_df3,by.x="Alyr_ID",by.y="Lyr_Gene",all.y=T)

Mias_arenosa_overlap_strict_df3<-Mias_arenosa_overlap_strict_df2[,-c(19:40)]
Mias_arenosa_overlap_strict_df4<-merge(Desc_HM,Mias_arenosa_overlap_strict_df3,by.x="Alyr_ID",by.y="Lyr_Gene",all.y=T)

Piek_arenosa_overlap_strict_df3<-Piek_arenosa_overlap_strict_df2[,-c(19:40)]
Piek_arenosa_overlap_strict_df4<-merge(Desc_HM,Piek_arenosa_overlap_strict_df3,by.x="Alyr_ID",by.y="Lyr_Gene",all.y=T)

write.xlsx(Kato_arenosa_overlap_strict_df4,"CNVs_overlap_strict_KatoKowa_arenosa_Lumpy_newortho.xlsx",sheetName="Kato_Kowa_arenosa",row.names=F,overwrite=T)
write.xlsx(Buko_arenosa_overlap_strict_df4,"CNVs_overlap_strict_BukoKowa_arenosa_Lumpy_newortho.xlsx",sheetName="Buko_Kowa_arenosa",row.names=F,overwrite=T)
write.xlsx(Mias_arenosa_overlap_strict_df4,"CNVs_overlap_strict_MiasKowa_arenosa_Lumpy_newortho.xlsx",sheetName="Mias_Kowa_arenosa",row.names=F,overwrite=T)
write.xlsx(Piek_arenosa_overlap_strict_df4,"CNVs_overlap_strict_PiekKowa_arenosa_Lumpy_newortho.xlsx",sheetName="Piek_Kowa_arenosa",row.names=F,overwrite=T)


write.xlsx(Kato_arenosa_overlap_strict_df4[!duplicated(Kato_arenosa_overlap_strict_df4$Alyr_ID),],"CNVs_overlap_strict_KatoKowa_arenosa_Lumpy_dedup_newortho.xlsx",sheetName="Kato_Kowa_arenosa",row.names=F,overwrite=T)
write.xlsx(Buko_arenosa_overlap_strict_df4[!duplicated(Buko_arenosa_overlap_strict_df4$Alyr_ID),],"CNVs_overlap_strict_BukoKowa_arenosa_Lumpy_dedup_newortho.xlsx",sheetName="Buko_Kowa_arenosa",row.names=F,overwrite=T)
write.xlsx(Mias_arenosa_overlap_strict_df4[!duplicated(Mias_arenosa_overlap_strict_df4$Alyr_ID),],"CNVs_overlap_strict_MiasKowa_arenosa_Lumpy_dedup_newortho.xlsx",sheetName="Mias_Kowa_arenosa",row.names=F,overwrite=T)
write.xlsx(Piek_arenosa_overlap_strict_df4[!duplicated(Piek_arenosa_overlap_strict_df4$Alyr_ID),],"CNVs_overlap_strict_PiekKowa_arenosa_Lumpy_dedup_newortho.xlsx",sheetName="Piek_Kowa_arenosa",row.names=F,overwrite=T)

nrow(Kato_arenosa_overlap_strict_df2)
#8839
nrow(Buko_arenosa_overlap_strict_df2)
#3628
nrow(Mias_arenosa_overlap_strict_df2)
#5478
nrow(Piek_arenosa_overlap_strict_df2)
#8212

length(unique(Kato_arenosa_overlap_strict_df2$Lyr_Gene))
#468
length(unique(Buko_arenosa_overlap_strict_df2$Lyr_Gene))
#355
length(unique(Mias_arenosa_overlap_strict_df2$Lyr_Gene))
#425
length(unique(Piek_arenosa_overlap_strict_df2$Lyr_Gene))
#469

#############
#Convergence#
#############
require(openxlsx)
Mias<-read.xlsx("CNVs_overlap_strict_MiasKowa_arenosa_Lumpy_dedup_newortho.xlsx",1)
Buko<-read.xlsx("CNVs_overlap_strict_BukoKowa_arenosa_Lumpy_dedup_newortho.xlsx",1)
Kato<-read.xlsx("CNVs_overlap_strict_KatoKowa_arenosa_Lumpy_dedup_newortho.xlsx",1)
Piek<-read.xlsx("CNVs_overlap_strict_PiekKowa_arenosa_Lumpy_dedup_newortho.xlsx",1)

MP<-merge(Mias,Piek,by="Alyr_ID")
nrow(MP)
#227

MPB<-merge(MP,Buko,by="Alyr_ID")
nrow(MPB)
#129

All<-merge(MPB,Kato,by="Alyr_ID")
nrow(All)
#107

write.xlsx(MP,"CNV_MiasPiek.xlsx",row.names=F)
write.xlsx(MPB,"CNV_MiasPiekBuko.xlsx",row.names=F)
write.xlsx(All,"CNV_All.xlsx",row.names=F)

MPnotKato<-MP[!MP$Alyr_ID%in%Kato$Alyr_ID,]
nrow(MPnotKato)
#66
write.xlsx(MPnotKato,"CNV_MiasPieknotKato.xlsx",row.names=F)




