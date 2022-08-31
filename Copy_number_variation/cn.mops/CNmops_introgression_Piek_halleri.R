#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesPiekhalleri <- scan("Piek_halleri_full.list",character(),quote="")
BAMFilesKowahalleri <- scan("Kowa_halleri_full.list",character(),quote="")

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
test<-scanBamHeader(BAMFilesPiekhalleri)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Piekhalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesPiekhalleri,sampleNames=c("Piek_h01","Piek_h03","Piek_h04","Piek_h05comb","Piek_h06comb","Piek_h07","Piek_h10","Piek_h11","Piek_h12","Piek_h13","Piek_h14_2","Piek_h15","Piek_h2"),refSeqName=sort(scaffolds),WL=500) 
Kowahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesKowahalleri,sampleNames=c("Kowa002h13","Kowa002h14","Kowa002h17","Kowa002h18","Kowa002h21","Kowa_h02comb","Kowa_h04comb","Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43"),refSeqName=sort(scaffolds),WL=500) 
PKha <- referencecn.mops(cases=Piekhalleri,controls=Kowahalleri, minWidth = 2, minReadCount=6)
resCNMOPSPKha <- calcIntegerCopyNumbers(PKha)

KPha <- referencecn.mops(cases=Kowahalleri,controls=Piekhalleri, minWidth = 2, minReadCount=6)
resCNMOPSKPha <- calcIntegerCopyNumbers(KPha)


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
regionsPKha<-cnvr(resCNMOPSPKha)
singlePKha<-cnvs(resCNMOPSPKha)
testPKha<-findOverlaps(regionsPKha,singlePKha)
test2PKha<- DataFrame(splitAsList(singlePKha$sampleName[subjectHits(testPKha)], queryHits(testPKha)),splitAsList(singlePKha$median[subjectHits(testPKha)], queryHits(testPKha)),splitAsList(singlePKha$mean[subjectHits(testPKha)], queryHits(testPKha)),splitAsList(singlePKha$CN[subjectHits(testPKha)], queryHits(testPKha)))
colnames(test2PKha)<-c("sampleNames","mean","median","CN")
test7PKha<-apply(as.data.frame(test2PKha),1,dupl_cov)
mcols(regionsPKha)<-as.data.frame(as.matrix(test7PKha))

dataPKha<-as.data.frame(regionsPKha)

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

adddataPKha<-do.call("rbind",apply(dataPKha,1,data_change))
dataPKha<-dataPKha[,-6]
frequPKha<-cbind(dataPKha,adddataPKha)

names(frequPKha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3PKha<-frequPKha[count.fields(textConnection(frequPKha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3PKha,"CnvregionswithmeanhalleriPiekKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigPKha<-read.xlsx("CnvregionswithmeanhalleriPiekKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigPKha[,2]<-as.numeric(as.character(frequbigPKha[,2]))
frequbigPKha[,3]<-as.numeric(as.character(frequbigPKha[,3]))

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
frequbigPKha_GRange<-GRanges(seqnames=tolower(frequbigPKha$Scaffold),ranges=IRanges(start=frequbigPKha$start,end=frequbigPKha$end))
values(frequbigPKha_GRange)<-frequbigPKha[,6:11]

testPKha<-findOverlaps(frequbigPKha_GRange,lyrGenes)
test2PKha<-pintersect(frequbigPKha_GRange[queryHits(testPKha)],lyrGenes[subjectHits(testPKha)])
test3PKha<-as.data.frame(test2PKha)
frequbigPKha_lyrGenes=mergeByOverlaps(frequbigPKha_GRange,lyrGenes,type=c("any"))
frequbigPKha_lyrGenes_df=data.frame(as.data.frame(frequbigPKha_lyrGenes$frequbigPKha_GRange),as.data.frame(frequbigPKha_lyrGenes$lyrGenes),test3PKha$width)
colnames(frequbigPKha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigPKha_lyrGenes_df<-frequbigPKha_lyrGenes_df[(frequbigPKha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigPKha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigPKha_lyrGenes_df,"PKha_temp.xlsx",col.names=TRUE,row.names=FALSE)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigPKha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigPKha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigPKha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup=merge(frequbigPKha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PKha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PKha.xlsx",overwrite=T)

Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigPKha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriPiekKowahallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Piek halleri#
#################################################

regionsKPha<-cnvr(resCNMOPSKPha)
singleKPha<-cnvs(resCNMOPSKPha)
testKPha<-findOverlaps(regionsKPha,singleKPha)
test2KPha<- DataFrame(splitAsList(singleKPha$sampleName[subjectHits(testKPha)], queryHits(testKPha)),splitAsList(singleKPha$median[subjectHits(testKPha)], queryHits(testKPha)),splitAsList(singleKPha$mean[subjectHits(testKPha)], queryHits(testKPha)),splitAsList(singleKPha$CN[subjectHits(testKPha)], queryHits(testKPha)))
colnames(test2KPha)<-c("sampleNames","mean","median","CN")
test7KPha<-apply(as.data.frame(test2KPha),1,dupl_cov)
mcols(regionsKPha)<-as.data.frame(as.matrix(test7KPha))

dataKPha<-as.data.frame(regionsKPha)

adddataKPha<-do.call("rbind",apply(dataKPha,1,data_change))
dataKPha<-dataKPha[,-6]
frequKPha<-cbind(dataKPha,adddataKPha)

names(frequKPha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KPha<-frequKPha[count.fields(textConnection(frequKPha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KPha,"CnvregionswithmeanhalleriKowaPiekhalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKPha<-read.xlsx("CnvregionswithmeanhalleriKowaPiekhalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigKPha[,2]<-as.numeric(as.character(frequbigKPha[,2]))
frequbigKPha[,3]<-as.numeric(as.character(frequbigKPha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKPha_GRange<-GRanges(seqnames=tolower(frequbigKPha$Scaffold),ranges=IRanges(start=frequbigKPha$start,end=frequbigKPha$end))
values(frequbigKPha_GRange)<-frequbigKPha[,6:11]

testKPha<-findOverlaps(frequbigKPha_GRange,lyrGenes)
test2KPha<-pintersect(frequbigKPha_GRange[queryHits(testKPha)],lyrGenes[subjectHits(testKPha)])
test3KPha<-as.data.frame(test2KPha)
frequbigKPha_lyrGenes=mergeByOverlaps(frequbigKPha_GRange,lyrGenes,type=c("any"))
frequbigKPha_lyrGenes_df=data.frame(as.data.frame(frequbigKPha_lyrGenes$frequbigKPha_GRange),as.data.frame(frequbigKPha_lyrGenes$lyrGenes),test3KPha$width)
colnames(frequbigKPha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKPha_lyrGenes_df<-frequbigKPha_lyrGenes_df[(frequbigKPha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKPha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKPha_lyrGenes_df,"KPha_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigKPha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKPha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKPha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup=merge(frequbigKPha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KPha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KPha.xlsx",overwrite=T)

Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKPha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriKowaPiekhallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

PKha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PKha.xlsx",1)
KPha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KPha.xlsx",1)

PKha_unique<-PKha[!duplicated(PKha$Lyr_Gene,PKha$Copy_start,PKha$Copy_end),]
KPha_unique<-KPha[!duplicated(KPha$Lyr_Gene,KPha$Copy_start,KPha$Copy_end),]

PKha_unique$Lyr_Gene<-as.character(PKha_unique$Lyr_Gene)
KPha_unique$Lyr_Gene<-as.character(KPha_unique$Lyr_Gene)
PKha_unique$CN_class<-as.character(PKha_unique$CN_class)
KPha_unique$CN_class<-as.character(KPha_unique$CN_class)

PKhaonly<-PKha_unique[!((PKha_unique$Lyr_Gene%in%KPha_unique$Lyr_Gene)&(PKha_unique$CN_class==KPha_unique$CN_class[match(PKha_unique$Lyr_Gene,KPha_unique$Lyr_Gene)])),]
KPhaonly<-KPha_unique[!((KPha_unique$Lyr_Gene%in%PKha_unique$Lyr_Gene)&(KPha_unique$CN_class==PKha_unique$CN_class[match(KPha_unique$Lyr_Gene,PKha_unique$Lyr_Gene)])),]

write.table(PKhaonly,"PKha_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(KPhaonly,"KPha_cnvs_only_strict.table",row.names=F,sep="\t")











