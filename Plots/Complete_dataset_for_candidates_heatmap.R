########################
#Load data#
########################

################
#FST#
################
setwd("F:/Introgression/Genome_scans/")
filelist = list.files(pattern = glob2rx("Fst_Mias_Kowa_arenosa_scaffold*_25SNPwins.table")) 
require(gtools)
filelist <- mixedsort(filelist)
MK_FST = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t")))
MK_FST<-MK_FST[MK_FST$End_position-MK_FST$Start_position<100000,]

filelist = list.files(pattern = glob2rx("Fst_Piek_Kowa_arenosa_scaffold*_25SNPwins.table")) 
filelist <- mixedsort(filelist)
PK_FST = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t")))
PK_FST<-PK_FST[PK_FST$End_position-PK_FST$Start_position<100000,]

filelist = list.files(pattern = glob2rx("Fst_Buko_Kowa_arenosa_scaffold*_25SNPwins.table")) 
filelist <- mixedsort(filelist)
BK_FST = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t")))
BK_FST<-BK_FST[BK_FST$End_position-BK_FST$Start_position<100000,]

filelist = list.files(pattern = glob2rx("Fst_Kato_Kowa_arenosa_scaffold*_25SNPwins.table")) 
filelist <- mixedsort(filelist)
KK_FST = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t")))
KK_FST<-KK_FST[KK_FST$End_position-KK_FST$Start_position<100000,]

filelist = list.files(pattern = glob2rx("Fst_Zapa_Kowa_arenosa_scaffold*_25SNPwins.table")) 
filelist <- mixedsort(filelist)
ZK_FST = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t")))
ZK_FST<-ZK_FST[ZK_FST$End_position-ZK_FST$Start_position<100000,]

filelist = list.files(pattern = glob2rx("Fst_Mias_arenosa_allhalleri_scaffold*_25SNPwins.table")) 
filelist <- mixedsort(filelist)
Mha_FST = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t")))
Mha_FST<-Mha_FST[Mha_FST$End_position-Mha_FST$Start_position<100000,]


################
#Twisst#
################

setwd("F:/Introgression/Twisst/Twisst_MiashalleriKowaarenosa_w100/")
chromNames = paste0("scaffold_",1:8)
weights_files4 = paste0("MiashalleriKowaarenosa.weights.w100.m1.",chromNames,".csv") 
window_data_files = paste0("MiashalleriKowaarenosalyr.Vcf.scf",chromNames,".w100m1.LDphase.phyml_bionj.data.tsv")
#make list of weight and window data for each chromsome, throwing out windows with NA values
weights_by_chrom4 = list()
window_data_by_chrom4 = list()
for (i in 1:length(weights_files4)){
    weights <- read.table(weights_files4[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom4[[i]] <- weights[good_rows,]
    window_data_by_chrom4[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }
weights_all4=data.frame()
for(i in 1:8){
weights_all4=rbind(weights_all4,weights_by_chrom4[[i]])}
mean_weights4 <- apply(weights_all4, 2, mean)
MK_twisst=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom4[[i]],weights_by_chrom4[[i]])
    MK_twisst=rbind(MK_twisst,df)
}


setwd("F:/Introgression/Twisst/Twisst_PiekhalleriKowaar_w100")
weights_files4 = paste0("PiekhalleriKowaarenosa.weights.w100.m1.",chromNames,".csv") 
window_data_files = paste0("PiekhalleriKowaarenosalyr.Vcf.scf",chromNames,".w100m1.LDphase.phyml_bionj.data.tsv")
weights_by_chrom4 = list()
window_data_by_chrom4 = list()
for (i in 1:length(weights_files4)){
    weights <- read.table(weights_files4[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom4[[i]] <- weights[good_rows,]
    window_data_by_chrom4[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }
weights_all4=data.frame()
for(i in 1:8){
weights_all4=rbind(weights_all4,weights_by_chrom4[[i]])}
mean_weights4 <- apply(weights_all4, 2, mean)
PK_twisst=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom4[[i]],weights_by_chrom4[[i]])
    PK_twisst=rbind(PK_twisst,df)
}


setwd("F:/Introgression/Twisst/Twisst_BukohalleriKowaar_w100")
weights_files4 = paste0("BukohalleriKowaarenosa.weights.w100.m1.",chromNames,".csv") 
window_data_files = paste0("BukohalleriKowaarenosalyr.Vcf.scf",chromNames,".w100m1.LDphase.phyml_bionj.data.tsv")
weights_by_chrom4 = list()
window_data_by_chrom4 = list()
for (i in 1:length(weights_files4)){
    weights <- read.table(weights_files4[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom4[[i]] <- weights[good_rows,]
    window_data_by_chrom4[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }
weights_all4=data.frame()
for(i in 1:8){
weights_all4=rbind(weights_all4,weights_by_chrom4[[i]])}
mean_weights4 <- apply(weights_all4, 2, mean)
BK_twisst=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom4[[i]],weights_by_chrom4[[i]])
    BK_twisst=rbind(BK_twisst,df)
}


setwd("F:/Introgression/Twisst/Twisst_KatohalleriKowaar_w100")
weights_files4 = paste0("KatohalleriKowaarenosa.weights.w100.m1.",chromNames,".csv") 
window_data_files = paste0("KatohalleriKowaarenosalyr.Vcf.scf",chromNames,".w100m1.LDphase.phyml_bionj.data.tsv")
weights_by_chrom4 = list()
window_data_by_chrom4 = list()
for (i in 1:length(weights_files4)){
    weights <- read.table(weights_files4[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom4[[i]] <- weights[good_rows,]
    window_data_by_chrom4[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }
weights_all4=data.frame()
for(i in 1:8){
weights_all4=rbind(weights_all4,weights_by_chrom4[[i]])}
mean_weights4 <- apply(weights_all4, 2, mean)
KK_twisst=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom4[[i]],weights_by_chrom4[[i]])
    KK_twisst=rbind(KK_twisst,df)
}


setwd("F:/Introgression/Twisst/Twisst_ZapahalleriKowaar_w100")
weights_files4 = paste0("ZapahalleriKowaarenosa.weights.w100.m1.",chromNames,".csv") 
window_data_files = paste0("ZapahalleriKowaarenosalyr.Vcf.scf",chromNames,".w100m1.LDphase.phyml_bionj.data.tsv")
weights_by_chrom4 = list()
window_data_by_chrom4 = list()
for (i in 1:length(weights_files4)){
    weights <- read.table(weights_files4[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom4[[i]] <- weights[good_rows,]
    window_data_by_chrom4[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }
weights_all4=data.frame()
for(i in 1:8){
weights_all4=rbind(weights_all4,weights_by_chrom4[[i]])}
mean_weights4 <- apply(weights_all4, 2, mean)
ZK_twisst=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom4[[i]],weights_by_chrom4[[i]])
    ZK_twisst=rbind(ZK_twisst,df)
}




#######################
#CNV#
#######################
setwd("F:/Introgression/Cnmops/Kowa/")
require(openxlsx)

MK_indCall<-read.xlsx("IndividualCall_data_cnmops_MKar.xlsx",1)
MK_normcov<-read.xlsx("Normalized_data_cnmops_MKar.xlsx",1)
mean_MK_indCall<-apply(MK_indCall,1,mean)
mean_MK_normcov<-apply(MK_normcov,1,mean)
MK<-data.frame(mean_MK_indCall,mean_MK_normcov)
Pos<-read.xlsx("Pos__data_cnmops_MKar.xlsx",1)
require(tidyr)
Pos<-separate(Pos,1,sep="_",into=c("NA","Scaff","Start","End"),remove=T)
Pos$Scaffold<-paste(Pos[,1],Pos[,2],sep="_")
Pos<-Pos[,c(5,3,4)]
Pos[,2]<-as.numeric(Pos[,2])
Pos[,3]<-as.numeric(Pos[,3])

MK_CNV<-data.frame(Pos,MK)

PK_indCall<-read.xlsx("IndividualCall_data_cnmops_PKar.xlsx",1)
PK_normcov<-read.xlsx("Normalized_data_cnmops_PKar.xlsx",1)
mean_PK_indCall<-apply(PK_indCall,1,mean)
mean_PK_normcov<-apply(PK_normcov,1,mean)
PK<-data.frame(mean_PK_indCall,mean_PK_normcov)
PK_CNV<-data.frame(Pos,PK)

BK_indCall<-read.xlsx("IndividualCall_data_cnmops_BKar.xlsx",1)
BK_normcov<-read.xlsx("Normalized_data_cnmops_BKar.xlsx",1)
mean_BK_indCall<-apply(BK_indCall,1,mean)
mean_BK_normcov<-apply(BK_normcov,1,mean)
BK<-data.frame(mean_BK_indCall,mean_BK_normcov)
BK_CNV<-data.frame(Pos,BK)

KK_indCall<-read.xlsx("IndividualCall_data_cnmops_KKar.xlsx",1)
KK_normcov<-read.xlsx("Normalized_data_cnmops_KKar.xlsx",1)
mean_KK_indCall<-apply(KK_indCall,1,mean)
mean_KK_normcov<-apply(KK_normcov,1,mean)
KK<-data.frame(mean_KK_indCall,mean_KK_normcov)
KK_CNV<-data.frame(Pos,KK)

ZK_indCall<-read.xlsx("IndividualCall_data_cnmops_ZKar.xlsx",1)
ZK_normcov<-read.xlsx("Normalized_data_cnmops_ZKar.xlsx",1)
mean_ZK_indCall<-apply(ZK_indCall,1,mean)
mean_ZK_normcov<-apply(ZK_normcov,1,mean)
ZK<-data.frame(mean_ZK_indCall,mean_ZK_normcov)
ZK_CNV<-data.frame(Pos,ZK)


#######################
#RNASeq#
#######################
setwd("F:/Introgression/Introrna/DESeq2/")
RNASeq<-read.table("Countdata_medofratios.table",sep="\t",header=T)


######################################
#merge to genes
######################################
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

## Load GFF files and add promotor region to genes

lyr_txdm = makeTxDbFromGFF(file ="F:/Cu_project/Lyrata/Alyrata_384_v2.1.gene.gff3", format = c("gff3"), organism = "Arabidopsis lyrata")

#lyr_txdm
lyrGenes=genes(lyr_txdm)
names(lyrGenes)<-NULL

#Define a function that allows the addition of promotor regions
extend <- function(x, upstream=0, downstream=0) #from https://support.bioconductor.org/p/78652/
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

# Use the function to extend promotor region 2kb
  # Upstream=2000 adds 2kb upstream from transcription start site
  # Downstream=0 adds 0bp downstream from the gene.
LyrGenes_plusPromot = extend(lyrGenes,upstream=2000, downstream=0)
LyrGenes_plusPromot<-LyrGenes_plusPromot[-which(lyrGenes@ elementMetadata@ listData$ gene_id=="AL1G54170"),]

MK_FST_GRange<-GRanges(seqnames=tolower(MK_FST$Scaffold),ranges=IRanges(start=MK_FST$Start_position,end=MK_FST$End_position))
values(MK_FST_GRange)<-MK_FST[,4]
MK_FST_genes=mergeByOverlaps(MK_FST_GRange,LyrGenes_plusPromot,type=c("any"))
MK_FST_genes_df=data.frame(as.data.frame(MK_FST_genes$MK_FST_GRange),as.data.frame(MK_FST_genes$LyrGenes_plusPromot))
MK_FST_genes_df=MK_FST_genes_df[order(MK_FST_genes_df$X,decreasing=T),c(1:3,6,12)]
colnames(MK_FST_genes_df)<-c("Scaffold","Window_start_FST_MK","Window_end_FST_MK","FST_MK","Alyr_ID")
MK_FST_genes_df=MK_FST_genes_df[!duplicated(MK_FST_genes_df$Alyr_ID),]

PK_FST_GRange<-GRanges(seqnames=tolower(PK_FST$Scaffold),ranges=IRanges(start=PK_FST$Start_position,end=PK_FST$End_position))
values(PK_FST_GRange)<-PK_FST[,4]
PK_FST_genes=mergeByOverlaps(PK_FST_GRange,LyrGenes_plusPromot,type=c("any"))
PK_FST_genes_df=data.frame(as.data.frame(PK_FST_genes$PK_FST_GRange),as.data.frame(PK_FST_genes$LyrGenes_plusPromot))
PK_FST_genes_df=PK_FST_genes_df[order(PK_FST_genes_df$X,decreasing=T),c(1:3,6,12)]
colnames(PK_FST_genes_df)<-c("Scaffold","Window_start_FST_PK","Window_end_FST_PK","FST_PK","Alyr_ID")
PK_FST_genes_df=PK_FST_genes_df[!duplicated(PK_FST_genes_df$Alyr_ID),]

BK_FST_GRange<-GRanges(seqnames=tolower(BK_FST$Scaffold),ranges=IRanges(start=BK_FST$Start_position,end=BK_FST$End_position))
values(BK_FST_GRange)<-BK_FST[,4]
BK_FST_genes=mergeByOverlaps(BK_FST_GRange,LyrGenes_plusPromot,type=c("any"))
BK_FST_genes_df=data.frame(as.data.frame(BK_FST_genes$BK_FST_GRange),as.data.frame(BK_FST_genes$LyrGenes_plusPromot))
BK_FST_genes_df=BK_FST_genes_df[order(BK_FST_genes_df$X,decreasing=T),c(1:3,6,12)]
colnames(BK_FST_genes_df)<-c("Scaffold","Window_start_FST_BK","Window_end_FST_BK","FST_BK","Alyr_ID")
BK_FST_genes_df=BK_FST_genes_df[!duplicated(BK_FST_genes_df$Alyr_ID),]

KK_FST_GRange<-GRanges(seqnames=tolower(KK_FST$Scaffold),ranges=IRanges(start=KK_FST$Start_position,end=KK_FST$End_position))
values(KK_FST_GRange)<-KK_FST[,4]
KK_FST_genes=mergeByOverlaps(KK_FST_GRange,LyrGenes_plusPromot,type=c("any"))
KK_FST_genes_df=data.frame(as.data.frame(KK_FST_genes$KK_FST_GRange),as.data.frame(KK_FST_genes$LyrGenes_plusPromot))
KK_FST_genes_df=KK_FST_genes_df[order(KK_FST_genes_df$X,decreasing=T),c(1:3,6,12)]
colnames(KK_FST_genes_df)<-c("Scaffold","Window_start_FST_KK","Window_end_FST_KK","FST_KK","Alyr_ID")
KK_FST_genes_df=KK_FST_genes_df[!duplicated(KK_FST_genes_df$Alyr_ID),]

ZK_FST_GRange<-GRanges(seqnames=tolower(ZK_FST$Scaffold),ranges=IRanges(start=ZK_FST$Start_position,end=ZK_FST$End_position))
values(ZK_FST_GRange)<-ZK_FST[,4]
ZK_FST_genes=mergeByOverlaps(ZK_FST_GRange,LyrGenes_plusPromot,type=c("any"))
ZK_FST_genes_df=data.frame(as.data.frame(ZK_FST_genes$ZK_FST_GRange),as.data.frame(ZK_FST_genes$LyrGenes_plusPromot))
ZK_FST_genes_df=ZK_FST_genes_df[order(ZK_FST_genes_df$X,decreasing=T),c(1:3,6,12)]
colnames(ZK_FST_genes_df)<-c("Scaffold","Window_start_FST_ZK","Window_end_FST_ZK","FST_ZK","Alyr_ID")
ZK_FST_genes_df=ZK_FST_genes_df[!duplicated(ZK_FST_genes_df$Alyr_ID),]


MK_twisst_GRange<-GRanges(seqnames=tolower(MK_twisst$scaffold),ranges=IRanges(start=MK_twisst$start,end=MK_twisst$end))
values(MK_twisst_GRange)<-MK_twisst[,c(7,9)]
MK_twisst_genes=mergeByOverlaps(MK_twisst_GRange,LyrGenes_plusPromot,type=c("any"))
MK_twisst_genes_df=data.frame(as.data.frame(MK_twisst_genes$MK_twisst_GRange),as.data.frame(MK_twisst_genes$LyrGenes_plusPromot))
MK_twisst_genes_df=MK_twisst_genes_df[order(MK_twisst_genes_df$topo3,decreasing=T),c(1:3,6:7,13)]
colnames(MK_twisst_genes_df)<-c("Scaffold","Window_start_Twisst_MK","Window_end_Twisst_MK","Twisst_topo1_MK","Twisst_topo3_MK","Alyr_ID")
MK_twisst_genes_df=MK_twisst_genes_df[!duplicated(MK_twisst_genes_df$Alyr_ID),]

PK_twisst_GRange<-GRanges(seqnames=tolower(PK_twisst$scaffold),ranges=IRanges(start=PK_twisst$start,end=PK_twisst$end))
values(PK_twisst_GRange)<-PK_twisst[,c(7,9)]
PK_twisst_genes=mergeByOverlaps(PK_twisst_GRange,LyrGenes_plusPromot,type=c("any"))
PK_twisst_genes_df=data.frame(as.data.frame(PK_twisst_genes$PK_twisst_GRange),as.data.frame(PK_twisst_genes$LyrGenes_plusPromot))
PK_twisst_genes_df=PK_twisst_genes_df[order(PK_twisst_genes_df$topo3,decreasing=T),c(1:3,6:7,13)]
colnames(PK_twisst_genes_df)<-c("Scaffold","Window_start_Twisst_PK","Window_end_Twisst_PK","Twisst_topo1_PK","Twisst_topo3_PK","Alyr_ID")
PK_twisst_genes_df=PK_twisst_genes_df[!duplicated(PK_twisst_genes_df$Alyr_ID),]

BK_twisst_GRange<-GRanges(seqnames=tolower(BK_twisst$scaffold),ranges=IRanges(start=BK_twisst$start,end=BK_twisst$end))
values(BK_twisst_GRange)<-BK_twisst[,c(7,9)]
BK_twisst_genes=mergeByOverlaps(BK_twisst_GRange,LyrGenes_plusPromot,type=c("any"))
BK_twisst_genes_df=data.frame(as.data.frame(BK_twisst_genes$BK_twisst_GRange),as.data.frame(BK_twisst_genes$LyrGenes_plusPromot))
BK_twisst_genes_df=BK_twisst_genes_df[order(BK_twisst_genes_df$topo3,decreasing=T),c(1:3,6:7,13)]
colnames(BK_twisst_genes_df)<-c("Scaffold","Window_start_Twisst_BK","Window_end_Twisst_BK","Twisst_topo1_BK","Twisst_topo3_BK","Alyr_ID")
BK_twisst_genes_df=BK_twisst_genes_df[!duplicated(BK_twisst_genes_df$Alyr_ID),]

KK_twisst_GRange<-GRanges(seqnames=tolower(KK_twisst$scaffold),ranges=IRanges(start=KK_twisst$start,end=KK_twisst$end))
values(KK_twisst_GRange)<-KK_twisst[,c(7,9)]
KK_twisst_genes=mergeByOverlaps(KK_twisst_GRange,LyrGenes_plusPromot,type=c("any"))
KK_twisst_genes_df=data.frame(as.data.frame(KK_twisst_genes$KK_twisst_GRange),as.data.frame(KK_twisst_genes$LyrGenes_plusPromot))
KK_twisst_genes_df=KK_twisst_genes_df[order(KK_twisst_genes_df$topo3,decreasing=T),c(1:3,6:7,13)]
colnames(KK_twisst_genes_df)<-c("Sccaffold","Window_start_Twisst_KK","Window_end_Twisst_KK","Twisst_topo1_KK","Twisst_topo3_KK","Alyr_ID")
KK_twisst_genes_df=KK_twisst_genes_df[!duplicated(KK_twisst_genes_df$Alyr_ID),]

ZK_twisst_GRange<-GRanges(seqnames=tolower(ZK_twisst$scaffold),ranges=IRanges(start=ZK_twisst$start,end=ZK_twisst$end))
values(ZK_twisst_GRange)<-ZK_twisst[,c(7,9)]
ZK_twisst_genes=mergeByOverlaps(ZK_twisst_GRange,LyrGenes_plusPromot,type=c("any"))
ZK_twisst_genes_df=data.frame(as.data.frame(ZK_twisst_genes$ZK_twisst_GRange),as.data.frame(ZK_twisst_genes$LyrGenes_plusPromot))
ZK_twisst_genes_df=ZK_twisst_genes_df[order(ZK_twisst_genes_df$topo3,decreasing=T),c(1:3,6:7,13)]
colnames(ZK_twisst_genes_df)<-c("Scaffold","Window_start_Twisst_ZK","Window_end_Twisst_ZK","Twisst_topo1_ZK","Twisst_topo3_ZK","Alyr_ID")
ZK_twisst_genes_df=ZK_twisst_genes_df[!duplicated(ZK_twisst_genes_df$Alyr_ID),]


MK_CNV_GRange<-GRanges(seqnames=tolower(MK_CNV$Scaffold),ranges=IRanges(start=MK_CNV$Start,end=MK_CNV$End))
values(MK_CNV_GRange)<-MK_CNV[,c(4:5)]
MK_CNV_genes=mergeByOverlaps(MK_CNV_GRange,LyrGenes_plusPromot,type=c("any"))
MK_CNV_genes_df=data.frame(as.data.frame(MK_CNV_genes$MK_CNV_GRange),as.data.frame(MK_CNV_genes$LyrGenes_plusPromot))
MK_CNV_genes_df=MK_CNV_genes_df[order(abs(MK_CNV_genes_df$mean_MK_indCall),decreasing=T),c(1:3,6:7,13)]
colnames(MK_CNV_genes_df)<-c("Scaffold","Window_start_CNV_MK","Window_end_CNV_MK","Cnmops_mean_indCall_MK","Cnmops_mean_normcov_MK","Alyr_ID")
MK_CNV_genes_df=MK_CNV_genes_df[!duplicated(MK_CNV_genes_df$Alyr_ID),]

PK_CNV_GRange<-GRanges(seqnames=tolower(PK_CNV$Scaffold),ranges=IRanges(start=PK_CNV$Start,end=PK_CNV$End))
values(PK_CNV_GRange)<-PK_CNV[,c(4:5)]
PK_CNV_genes=mergeByOverlaps(PK_CNV_GRange,LyrGenes_plusPromot,type=c("any"))
PK_CNV_genes_df=data.frame(as.data.frame(PK_CNV_genes$PK_CNV_GRange),as.data.frame(PK_CNV_genes$LyrGenes_plusPromot))
PK_CNV_genes_df=PK_CNV_genes_df[order(abs(PK_CNV_genes_df$mean_PK_indCall),decreasing=T),c(1:3,6:7,13)]
colnames(PK_CNV_genes_df)<-c("Scaffold","Window_start_CNV_PK","Window_end_CNV_PK","Cnmops_mean_indCall_PK","Cnmops_mean_normcov_PK","Alyr_ID")
PK_CNV_genes_df=PK_CNV_genes_df[!duplicated(PK_CNV_genes_df$Alyr_ID),]

BK_CNV_GRange<-GRanges(seqnames=tolower(BK_CNV$Scaffold),ranges=IRanges(start=BK_CNV$Start,end=BK_CNV$End))
values(BK_CNV_GRange)<-BK_CNV[,c(4:5)]
BK_CNV_genes=mergeByOverlaps(BK_CNV_GRange,LyrGenes_plusPromot,type=c("any"))
BK_CNV_genes_df=data.frame(as.data.frame(BK_CNV_genes$BK_CNV_GRange),as.data.frame(BK_CNV_genes$LyrGenes_plusPromot))
BK_CNV_genes_df=BK_CNV_genes_df[order(abs(BK_CNV_genes_df$mean_BK_indCall),decreasing=T),c(1:3,6:7,13)]
colnames(BK_CNV_genes_df)<-c("Scaffold","Window_start_CNV_BK","Window_end_CNV_BK","Cnmops_mean_indCall_BK","Cnmops_mean_normcov_BK","Alyr_ID")
BK_CNV_genes_df=BK_CNV_genes_df[!duplicated(BK_CNV_genes_df$Alyr_ID),]

KK_CNV_GRange<-GRanges(seqnames=tolower(KK_CNV$Scaffold),ranges=IRanges(start=KK_CNV$Start,end=KK_CNV$End))
values(KK_CNV_GRange)<-KK_CNV[,c(4:5)]
KK_CNV_genes=mergeByOverlaps(KK_CNV_GRange,LyrGenes_plusPromot,type=c("any"))
KK_CNV_genes_df=data.frame(as.data.frame(KK_CNV_genes$KK_CNV_GRange),as.data.frame(KK_CNV_genes$LyrGenes_plusPromot))
KK_CNV_genes_df=KK_CNV_genes_df[order(abs(KK_CNV_genes_df$mean_KK_indCall),decreasing=T),c(1:3,6:7,13)]
colnames(KK_CNV_genes_df)<-c("Scaffold","Window_start_CNV_KK","Window_end_CNV_KK","Cnmops_mean_indCall_KK","Cnmops_mean_normcov_KK","Alyr_ID")
KK_CNV_genes_df=KK_CNV_genes_df[!duplicated(KK_CNV_genes_df$Alyr_ID),]

ZK_CNV_GRange<-GRanges(seqnames=tolower(ZK_CNV$Scaffold),ranges=IRanges(start=ZK_CNV$Start,end=ZK_CNV$End))
values(ZK_CNV_GRange)<-ZK_CNV[,c(4:5)]
ZK_CNV_genes=mergeByOverlaps(ZK_CNV_GRange,LyrGenes_plusPromot,type=c("any"))
ZK_CNV_genes_df=data.frame(as.data.frame(ZK_CNV_genes$ZK_CNV_GRange),as.data.frame(ZK_CNV_genes$LyrGenes_plusPromot))
ZK_CNV_genes_df=ZK_CNV_genes_df[order(abs(ZK_CNV_genes_df$mean_ZK_indCall),decreasing=T),c(1:3,6:7,13)]
colnames(ZK_CNV_genes_df)<-c("Scaffold","Window_start_CNV_ZK","Window_end_CNV_ZK","Cnmops_mean_indCall_ZK","Cnmops_mean_normcov_ZK","Alyr_ID")
ZK_CNV_genes_df=ZK_CNV_genes_df[!duplicated(ZK_CNV_genes_df$Alyr_ID),]

MK_all<-merge(merge(MK_twisst_genes_df,MK_CNV_genes_df,by="Alyr_ID"),MK_FST_genes_df,by="Alyr_ID")
#29852
PK_all<-merge(merge(PK_twisst_genes_df,PK_CNV_genes_df,by="Alyr_ID"),PK_FST_genes_df,by="Alyr_ID")
#29943
BK_all<-merge(merge(BK_twisst_genes_df,BK_CNV_genes_df,by="Alyr_ID"),BK_FST_genes_df,by="Alyr_ID")
#29789
KK_all<-merge(merge(KK_twisst_genes_df,KK_CNV_genes_df,by="Alyr_ID"),KK_FST_genes_df,by="Alyr_ID")
#29847
ZK_all<-merge(merge(ZK_twisst_genes_df,ZK_CNV_genes_df,by="Alyr_ID"),ZK_FST_genes_df,by="Alyr_ID")
#29803

MK_all2<-MK_all[,c(1:6,8:11,13:15)]
colnames(MK_all2)<-c("Alyr_ID","Scaffold","Window_start_Twisst","Window_end_Twisst","Twisst_topo1","Twisst_topo3",
"Window_start_CNV","Window_end_CNV","Cnmops_mean_indCall","Cnmops_mean_normcov","Window_start_FST","Window_end_FST","FST")

PK_all2<-PK_all[,c(1:6,8:11,13:15)]
colnames(PK_all2)<-c("Alyr_ID","Scaffold","Window_start_Twisst","Window_end_Twisst","Twisst_topo1","Twisst_topo3",
"Window_start_CNV","Window_end_CNV","Cnmops_mean_indCall","Cnmops_mean_normcov","Window_start_FST","Window_end_FST","FST")

BK_all2<-BK_all[,c(1:6,8:11,13:15)]
colnames(BK_all2)<-c("Alyr_ID","Scaffold","Window_start_Twisst","Window_end_Twisst","Twisst_topo1","Twisst_topo3",
"Window_start_CNV","Window_end_CNV","Cnmops_mean_indCall","Cnmops_mean_normcov","Window_start_FST","Window_end_FST","FST")

KK_all2<-KK_all[,c(1:6,8:11,13:15)]
colnames(KK_all2)<-c("Alyr_ID","Scaffold","Window_start_Twisst","Window_end_Twisst","Twisst_topo1","Twisst_topo3",
"Window_start_CNV","Window_end_CNV","Cnmops_mean_indCall","Cnmops_mean_normcov","Window_start_FST","Window_end_FST","FST")

ZK_all2<-ZK_all[,c(1:6,8:11,13:15)]
colnames(ZK_all2)<-c("Alyr_ID","Scaffold","Window_start_Twisst","Window_end_Twisst","Twisst_topo1","Twisst_topo3",
"Window_start_CNV","Window_end_CNV","Cnmops_mean_indCall","Cnmops_mean_normcov","Window_start_FST","Window_end_FST","FST")

require(RColorBrewer)

require(vegan)

matrix_all<-data.frame(rbind(MK_all2[,c(1,5:6,9:10,13)],PK_all2[,c(1,5:6,9:10,13)],
BK_all2[,c(1,5:6,9:10,13)],KK_all2[,c(1,5:6,9:10,13)],
ZK_all2[,c(1,5:6,9:10,13)]))

matrix_all$Pop<-c(rep("Mias",nrow(MK_all2)),rep("Piek",nrow(PK_all2)),rep("Buko",nrow(BK_all2)),rep("Kato",nrow(KK_all2)),rep("Zapa",nrow(ZK_all2)))
summary(matrix_all)

matrix_all$Cnmops_mean_indCall<-ifelse(matrix_all$Cnmops_mean_indCall<(-2),-2,matrix_all$Cnmops_mean_indCall)
matrix_all$FST<-ifelse(matrix_all$FST<0,0,matrix_all$FST)
matrix_all$Cnmops_mean_indCall<-decostand(matrix_all$Cnmops_mean_indCall,"range")
matrix_all$Cnmops_mean_normcov<-decostand(matrix_all$Cnmops_mean_normcov,"range")
matrix_all$FST<-decostand(matrix_all$FST,"range")


#HMA4,MTP1,HMA3,ZIP6,NRAMP3,ZIP11
Candlist<-c("AL3G52820","AL4G46270","AL7G22080","AL4G24850","AL4G13010","AL1G64050")

matrix_cand4<-matrix(1:4,ncol=1)
for (j in c(6,4,3))
	{for (i in seq(6,1,-1))
		{Cand=Candlist[i]
		matrix_cand4<-as.matrix(cbind(matrix_cand4,rbind(matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Mias",j],matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Piek",j],matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Buko",j],matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Kato",j])))		
		}
	matrix_cand4<-as.matrix(cbind(matrix_cand4,rbind(median(matrix_all[matrix_all$Pop=="Mias",j]),median(matrix_all[matrix_all$Pop=="Piek",j]),median(matrix_all[matrix_all$Pop=="Buko",j]),median(matrix_all[matrix_all$Pop=="Kato",j])),NA))
	}
matrix_cand4<-matrix_cand4[,-1]
colnames(matrix_cand4)<-rep(rev(c("","Median","HMA4","MTP1","HMA3","ZIP6","NRAMP3","ZIP11")),3)
matrix_cand4<-matrix_cand4[,c(7,1:6,8,15,9:14,16,23,17:22,24)]
rownames(matrix_cand4)<-c("Mias","Piek","Buko","Kato")

pdf("All_cand_genes_heatmap_paper.pdf",width=10,height=20,paper="special",pointsize=20)
par(oma=c(1,1,1,1))
par(mar=c(1,1,1,1))
heatmap(t(matrix_cand4),col=colorRampPalette(brewer.pal(9, "Reds") )(255)
,ylab="",xlab="",Rowv=NA,Colv=NA,scale="none",margins=c(5,22),las=1)
#text(x=c(seq(0.1,0.9,length.out=3)-0.1),y=0.05,labels = c("Topo3 weighting","Mean CN call",expression(F[ST])),srt=45,cex=1.2,xpd=T,adj=1)

dev.off()

pdf("All_cand_genes_heatmap_paper_colourbar_legend.pdf",width=10,height=2,paper="special",pointsize=20)
par(mar=c(3,1,1,1))
image(as.matrix(1:100),col=colorRampPalette(brewer.pal(9, "Reds") )(255),ylab="",xlab="",yaxt="n",cex.axis=1.3,useRaster = TRUE)
dev.off()

#IRT2, ZIP8, FRO5, FER1, FRD3
Candlist2<-c("AL7G34360","AL8G14120","AL6G35310","AL6G10450","AL3G19430")

matrix_cand5<-matrix(1:4,ncol=1)
for (j in c(6,4,3))
	{for (i in seq(5,1,-1))
		{Cand=Candlist2[i]
		matrix_cand5<-as.matrix(cbind(matrix_cand5,rbind(matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Mias",j],matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Piek",j],matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Buko",j],matrix_all[matrix_all$Alyr_ID==Cand&matrix_all$Pop=="Kato",j])))		
		}
	matrix_cand5<-as.matrix(cbind(matrix_cand5,rbind(median(matrix_all[matrix_all$Pop=="Mias",j]),median(matrix_all[matrix_all$Pop=="Piek",j]),median(matrix_all[matrix_all$Pop=="Buko",j]),median(matrix_all[matrix_all$Pop=="Kato",j])),NA))
	}
matrix_cand5<-matrix_cand5[,-1]
colnames(matrix_cand5)<-rep(rev(c("","Median","IRT2","ZIP8","FRO5","FER1","FRD3")),3)
matrix_cand5<-matrix_cand5[,c(6,1:5,7,13,8:12,14,20,15:19,21)]
rownames(matrix_cand5)<-c("Mias","Piek","Buko","Kato")

pdf("All_cand_genes_heatmap_paper_set2.pdf",width=10,height=20,paper="special",pointsize=20)
par(oma=c(1,1,1,1))
par(mar=c(1,1,1,1))
heatmap(t(matrix_cand5),col=colorRampPalette(brewer.pal(9, "Reds") )(255)
,ylab="",xlab="",Rowv=NA,Colv=NA,scale="none",margins=c(5,22),las=1)
#text(x=c(seq(0.1,0.9,length.out=3)-0.1),y=0.05,labels = c("Topo3 weighting","Mean CN call",expression(F[ST])),srt=45,cex=1.2,xpd=T,adj=1)

dev.off()

colorRampPalette(brewer.pal(9, "OrRd") )(255)
rev(colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255))
colorRampPalette(brewer.pal(9,"Greys"))(11)
colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)


