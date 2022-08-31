#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesBukohalleri <- scan("Buko_halleri_full.list",character(),quote="")
#BAMFilesBukoarenosa <- scan("Buko_arenosa_full.list",character(),quote="")
#BAMFilesKatohalleri <- scan("Kato_halleri_full.list",character(),quote="")
#BAMFilesKatoarenosa <- scan("Kato_arenosa_full.list",character(),quote="")
#BAMFilesMiashalleri <- scan("Mias_halleri_full.list",character(),quote="")
#BAMFilesMiasarenosa <- scan("Mias_arenosa_full.list",character(),quote="")
#BAMFilesPiekhalleri <- scan("Piek_halleri_full.list",character(),quote="")
#BAMFilesPiekarenosa <- scan("Piek_arenosa_full.list",character(),quote="")
BAMFilesKowahalleri <- scan("Kowa_halleri_full.list",character(),quote="")
#BAMFilesKowaarenosa <- scan("Kowa_arenosa_full.list",character(),quote="")
#BAMFilesZapahalleri <- scan("Zapa_halleri_full.list",character(),quote="")
#BAMFilesZapaarenosa <- scan("Zapa_arenosa_full.list",character(),quote="")


#Window length set to: 500
#C=m*RL/WL
#C: coverage
#m: average number of reads per bin ~ 50-100
#RL: Read length
#WL: window length

#WL=m*RL/C
#WL=100*125/19 ~658
#or WL=100*150/15~1000
#m should be between 50 and 100 https://support.bioconductor.org/p/75969/ 25.01.16

require(Rsamtools)
test<-scanBamHeader(BAMFilesBukohalleri)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Bukohalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesBukohalleri,sampleNames=c("Buko_h04","Buko_h12","Buko_h1_7","Buko_h21","Buko_h22","Buko_h23","Buko_h24","Buko_h25","Buko_h26","Buko_h27","Buko_h28","Buko_h29","Buko_h30","Buko_h31","Buko_h32","Buko_h33","Buko_h35","Buko_h36"),refSeqName=sort(scaffolds),WL=500) 
Kowahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesKowahalleri,sampleNames=c("Kowa002h13","Kowa002h14","Kowa002h17","Kowa002h18","Kowa002h21","Kowa_h02comb","Kowa_h04comb","Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43"),refSeqName=sort(scaffolds),WL=500) 
BKha <- referencecn.mops(cases=Bukohalleri,controls=Kowahalleri, minWidth = 2, minReadCount=6)
resCNMOPSBKha <- calcIntegerCopyNumbers(BKha)

#Bukoarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesBukoarenosa,sampleNames=c("Buko_h04","Buko_h12","Buko_h1_7","Buko_h21","Buko_h22","Buko_h23","Buko_h24","Buko_h25","Buko_h26","Buko_h27","Buko_h28","Buko_h29","Buko_h30","Buko_h31","Buko_h32","Buko_h33","Buko_h35","Buko_h36"),refSeqName=sort(scaffolds),WL=500) 
#Kowaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesKowahalleri,sampleNames=c("Kowa002h13","Kowa002h14","Kowa002h17","Kowa002h18","Kowa002h21","Kowa_h02comb","Kowa_h04comb","Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43"),refSeqName=sort(scaffolds),WL=500) 
#BKa <- referencecn.mops(cases=Bukoarenosa,controls=Kowaarenosa, minWidth = 2, minReadCount=6)
#resCNMOPSBKha <- calcIntegerCopyNumbers(BKa)


KBha <- referencecn.mops(cases=Kowahalleri,controls=Bukohalleri, minWidth = 2, minReadCount=6)
resCNMOPSKBha <- calcIntegerCopyNumbers(KBha)


#########################################################
#overlap between individual data and regions to get mean#
#########################################################
#after refrencecn.mops and CalcIntegerCopyNumbers
#########################################################
#overlap between individual data and regions to get mean#
#########################################################
dupl_cov<-function(arg1){
	arg2<-arg1
	arg2[[2]]<-round(arg2[[2]],digits=2)
	arg2[[3]]<-round(arg2[[3]],digits=2)
	if (any(duplicated(arg2[[1]]))){
	duplist<-c(as.character(arg2[[1]][duplicated(arg2[[1]])]))
	for (i in duplist){
		arg2[[2]][arg2[[1]]==i]<-mean(arg2[[2]][arg2[[1]]==i])
		arg2[[3]][arg2[[1]]==i]<-median(arg2[[3]][arg2[[1]]==i])
		arg2[[4]][arg2[[1]]==i]<-names(which.max(table(arg2[[4]][arg2[[1]]==i])))
		}
	arg2[[2]]<-arg2[[2]][!duplicated(arg2[[1]])]
	arg2[[3]]<-arg2[[3]][!duplicated(arg2[[1]])]
	arg2[[4]]<-arg2[[4]][!duplicated(arg2[[1]])]
	arg2[[1]]<-as.character(arg2[[1]][!duplicated(arg2[[1]])])
	}
	return(arg2)
}


#after refrencecn.mops and CalcIntegerCopyNumbers
regionsBKha<-cnvr(resCNMOPSBKha)
singleBKha<-cnvs(resCNMOPSBKha)
testBKha<-findOverlaps(regionsBKha,singleBKha)
test2BKha<- DataFrame(splitAsList(singleBKha$sampleName[subjectHits(testBKha)], queryHits(testBKha)),splitAsList(singleBKha$median[subjectHits(testBKha)], queryHits(testBKha)),splitAsList(singleBKha$mean[subjectHits(testBKha)], queryHits(testBKha)),splitAsList(singleBKha$CN[subjectHits(testBKha)], queryHits(testBKha)))
colnames(test2BKha)<-c("sampleNames","mean","median","CN")
test7BKha<-apply(as.data.frame(test2BKha),1,dupl_cov)
mcols(regionsBKha)<-as.data.frame(as.matrix(test7BKha))

dataBKha<-as.data.frame(regionsBKha)

data_change<-function(arg1){
	arg2<-data.frame(matrix(ncol=6,nrow=0))
	arg2[[1,1]]<-as.character(sapply(arg1[[6]][1],paste0,collapse=",")) 
	arg2[[1,2]]<-sapply(arg1[[6]][2],paste0,collapse=",") 
	arg2[[1,3]]<-sapply(arg1[[6]][3],paste0,collapse=",") 
	arg2[[1,4]]<-sapply(arg1[[6]][4],paste0,collapse=",") 
	arg2[[1,5]]<-names(which.max(table(arg1[[6]][4])))
	arg2[[1,6]]<-max(table(arg1[[6]][4]))/lengths(arg1[[6]][4])
	return(arg2)
}

adddataBKha<-do.call("rbind",apply(dataBKha,1,data_change))
dataBKha<-dataBKha[,-6]
frequBKha<-cbind(dataBKha,adddataBKha)

names(frequBKha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3BKha<-frequBKha[count.fields(textConnection(frequBKha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3BKha,"CnvregionswithmeanhalleriBukoKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigBKha<-read.xlsx("CnvregionswithmeanhalleriBukoKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigBKha[,2]<-as.numeric(as.character(frequbigBKha[,2]))
frequbigBKha[,3]<-as.numeric(as.character(frequbigBKha[,3]))

###################################################
#find overlaps between windows and A. lyrata genes#
###################################################
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

### Step 1: import genes from GFF files

lyr_txdm = makeTxDbFromGFF(file ="Alyrata_384_v2.1.gene.gff3", format = c("gff3"), organism = "Arabidopsis lyrata")

#lyr_txdm
lyrGenes=genes(lyr_txdm)
names(lyrGenes)<-NULL

### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigBKha_GRange<-GRanges(seqnames=tolower(frequbigBKha$Scaffold),ranges=IRanges(start=frequbigBKha$start,end=frequbigBKha$end))
values(frequbigBKha_GRange)<-frequbigBKha[,6:11]

testBKha<-findOverlaps(frequbigBKha_GRange,lyrGenes)
test2BKha<-pintersect(frequbigBKha_GRange[queryHits(testBKha)],lyrGenes[subjectHits(testBKha)])
test3BKha<-as.data.frame(test2BKha)
frequbigBKha_lyrGenes=mergeByOverlaps(frequbigBKha_GRange,lyrGenes,type=c("any"))
frequbigBKha_lyrGenes_df=data.frame(as.data.frame(frequbigBKha_lyrGenes$frequbigBKha_GRange),as.data.frame(frequbigBKha_lyrGenes$lyrGenes),test3BKha$width)
colnames(frequbigBKha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigBKha_lyrGenes_df<-frequbigBKha_lyrGenes_df[(frequbigBKha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigBKha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigBKha_lyrGenes_df,"BKha_temp.xlsx",col.names=TRUE,row.names=FALSE)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigBKha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigBKha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigBKha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup=merge(frequbigBKha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BKha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BKha.xlsx",overwrite=T)

Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigBKha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriBukoKowahallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Buko halleri#
#################################################

regionsKBha<-cnvr(resCNMOPSKBha)
singleKBha<-cnvs(resCNMOPSKBha)
testKBha<-findOverlaps(regionsKBha,singleKBha)
test2KBha<- DataFrame(splitAsList(singleKBha$sampleName[subjectHits(testKBha)], queryHits(testKBha)),splitAsList(singleKBha$median[subjectHits(testKBha)], queryHits(testKBha)),splitAsList(singleKBha$mean[subjectHits(testKBha)], queryHits(testKBha)),splitAsList(singleKBha$CN[subjectHits(testKBha)], queryHits(testKBha)))
colnames(test2KBha)<-c("sampleNames","mean","median","CN")
test7KBha<-apply(as.data.frame(test2KBha),1,dupl_cov)
mcols(regionsKBha)<-as.data.frame(as.matrix(test7KBha))

dataKBha<-as.data.frame(regionsKBha)

adddataKBha<-do.call("rbind",apply(dataKBha,1,data_change))
dataKBha<-dataKBha[,-6]
frequKBha<-cbind(dataKBha,adddataKBha)

names(frequKBha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KBha<-frequKBha[count.fields(textConnection(frequKBha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KBha,"CnvregionswithmeanhalleriKowaBukohalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKBha<-read.xlsx("CnvregionswithmeanhalleriKowaBukohalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigKBha[,2]<-as.numeric(as.character(frequbigKBha[,2]))
frequbigKBha[,3]<-as.numeric(as.character(frequbigKBha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKBha_GRange<-GRanges(seqnames=tolower(frequbigKBha$Scaffold),ranges=IRanges(start=frequbigKBha$start,end=frequbigKBha$end))
values(frequbigKBha_GRange)<-frequbigKBha[,6:11]

testKBha<-findOverlaps(frequbigKBha_GRange,lyrGenes)
test2KBha<-pintersect(frequbigKBha_GRange[queryHits(testKBha)],lyrGenes[subjectHits(testKBha)])
test3KBha<-as.data.frame(test2KBha)
frequbigKBha_lyrGenes=mergeByOverlaps(frequbigKBha_GRange,lyrGenes,type=c("any"))
frequbigKBha_lyrGenes_df=data.frame(as.data.frame(frequbigKBha_lyrGenes$frequbigKBha_GRange),as.data.frame(frequbigKBha_lyrGenes$lyrGenes),test3KBha$width)
colnames(frequbigKBha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKBha_lyrGenes_df<-frequbigKBha_lyrGenes_df[(frequbigKBha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKBha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKBha_lyrGenes_df,"KBha_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigKBha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKBha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKBha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup=merge(frequbigKBha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KBha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KBha.xlsx",overwrite=T)

Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKBha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriKowaBukohallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

BKha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BKha.xlsx",1)
KBha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KBha.xlsx",1)

BKha_unique<-BKha[!duplicated(BKha$Lyr_Gene,BKha$Copy_start,BKha$Copy_end),]
KBha_unique<-KBha[!duplicated(KBha$Lyr_Gene,KBha$Copy_start,KBha$Copy_end),]

BKha_unique$Lyr_Gene<-as.character(BKha_unique$Lyr_Gene)
KBha_unique$Lyr_Gene<-as.character(KBha_unique$Lyr_Gene)
BKha_unique$CN_class<-as.character(BKha_unique$CN_class)
KBha_unique$CN_class<-as.character(KBha_unique$CN_class)

BKhaonly<-BKha_unique[!((BKha_unique$Lyr_Gene%in%KBha_unique$Lyr_Gene)&(BKha_unique$CN_class==KBha_unique$CN_class[match(BKha_unique$Lyr_Gene,KBha_unique$Lyr_Gene)])),]
KBhaonly<-KBha_unique[!((KBha_unique$Lyr_Gene%in%BKha_unique$Lyr_Gene)&(KBha_unique$CN_class==BKha_unique$CN_class[match(KBha_unique$Lyr_Gene,BKha_unique$Lyr_Gene)])),]

write.table(BKhaonly,"BKha_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(KBhaonly,"KBha_cnvs_only_strict.table",row.names=F,sep="\t")











