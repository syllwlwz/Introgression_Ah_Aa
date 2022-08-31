#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesMiashalleri <- scan("Mias_halleri_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesMiashalleri)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Miashalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesMiashalleri,sampleNames=c("Mias002h19","Mias002h20","Mias004h02","Mias004h06","Mias004h07","Mias004h08","Mias009h03","Mias009h05","Mias_a54","Mias_h09","Mias_h41","Mias_h43","Mias_h44","Mias_h45","Mias_h46"),refSeqName=sort(scaffolds),WL=500) 
Kowahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesKowahalleri,sampleNames=c("Kowa002h13","Kowa002h14","Kowa002h17","Kowa002h18","Kowa002h21","Kowa_h02comb","Kowa_h04comb","Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43"),refSeqName=sort(scaffolds),WL=500) 
MKha <- referencecn.mops(cases=Miashalleri,controls=Kowahalleri, minWidth = 2, minReadCount=6)
resCNMOPSMKha <- calcIntegerCopyNumbers(MKha)

KMha <- referencecn.mops(cases=Kowahalleri,controls=Miashalleri, minWidth = 2, minReadCount=6)
resCNMOPSKMha <- calcIntegerCopyNumbers(KMha)


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
regionsMKha<-cnvr(resCNMOPSMKha)
singleMKha<-cnvs(resCNMOPSMKha)
testMKha<-findOverlaps(regionsMKha,singleMKha)
test2MKha<- DataFrame(splitAsList(singleMKha$sampleName[subjectHits(testMKha)], queryHits(testMKha)),splitAsList(singleMKha$median[subjectHits(testMKha)], queryHits(testMKha)),splitAsList(singleMKha$mean[subjectHits(testMKha)], queryHits(testMKha)),splitAsList(singleMKha$CN[subjectHits(testMKha)], queryHits(testMKha)))
colnames(test2MKha)<-c("sampleNames","mean","median","CN")
test7MKha<-apply(as.data.frame(test2MKha),1,dupl_cov)
mcols(regionsMKha)<-as.data.frame(as.matrix(test7MKha))

dataMKha<-as.data.frame(regionsMKha)

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

adddataMKha<-do.call("rbind",apply(dataMKha,1,data_change))
dataMKha<-dataMKha[,-6]
frequMKha<-cbind(dataMKha,adddataMKha)

names(frequMKha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3MKha<-frequMKha[count.fields(textConnection(frequMKha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3MKha,"CnvregionswithmeanhalleriMiasKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigMKha<-read.xlsx("CnvregionswithmeanhalleriMiasKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigMKha[,2]<-as.numeric(as.character(frequbigMKha[,2]))
frequbigMKha[,3]<-as.numeric(as.character(frequbigMKha[,3]))

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
frequbigMKha_GRange<-GRanges(seqnames=tolower(frequbigMKha$Scaffold),ranges=IRanges(start=frequbigMKha$start,end=frequbigMKha$end))
values(frequbigMKha_GRange)<-frequbigMKha[,6:11]

testMKha<-findOverlaps(frequbigMKha_GRange,lyrGenes)
test2MKha<-pintersect(frequbigMKha_GRange[queryHits(testMKha)],lyrGenes[subjectHits(testMKha)])
test3MKha<-as.data.frame(test2MKha)
frequbigMKha_lyrGenes=mergeByOverlaps(frequbigMKha_GRange,lyrGenes,type=c("any"))
frequbigMKha_lyrGenes_df=data.frame(as.data.frame(frequbigMKha_lyrGenes$frequbigMKha_GRange),as.data.frame(frequbigMKha_lyrGenes$lyrGenes),test3MKha$width)
colnames(frequbigMKha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigMKha_lyrGenes_df<-frequbigMKha_lyrGenes_df[(frequbigMKha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigMKha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigMKha_lyrGenes_df,"MKha_temp.xlsx",col.names=TRUE,row.names=FALSE)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigMKha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigMKha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigMKha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup=merge(frequbigMKha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MKha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MKha.xlsx",overwrite=T)

Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigMKha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriMiasKowahallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Mias halleri#
#################################################

regionsKMha<-cnvr(resCNMOPSKMha)
singleKMha<-cnvs(resCNMOPSKMha)
testKMha<-findOverlaps(regionsKMha,singleKMha)
test2KMha<- DataFrame(splitAsList(singleKMha$sampleName[subjectHits(testKMha)], queryHits(testKMha)),splitAsList(singleKMha$median[subjectHits(testKMha)], queryHits(testKMha)),splitAsList(singleKMha$mean[subjectHits(testKMha)], queryHits(testKMha)),splitAsList(singleKMha$CN[subjectHits(testKMha)], queryHits(testKMha)))
colnames(test2KMha)<-c("sampleNames","mean","median","CN")
test7KMha<-apply(as.data.frame(test2KMha),1,dupl_cov)
mcols(regionsKMha)<-as.data.frame(as.matrix(test7KMha))

dataKMha<-as.data.frame(regionsKMha)

adddataKMha<-do.call("rbind",apply(dataKMha,1,data_change))
dataKMha<-dataKMha[,-6]
frequKMha<-cbind(dataKMha,adddataKMha)

names(frequKMha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KMha<-frequKMha[count.fields(textConnection(frequKMha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KMha,"CnvregionswithmeanhalleriKowaMiashalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKMha<-read.xlsx("CnvregionswithmeanhalleriKowaMiashalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigKMha[,2]<-as.numeric(as.character(frequbigKMha[,2]))
frequbigKMha[,3]<-as.numeric(as.character(frequbigKMha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKMha_GRange<-GRanges(seqnames=tolower(frequbigKMha$Scaffold),ranges=IRanges(start=frequbigKMha$start,end=frequbigKMha$end))
values(frequbigKMha_GRange)<-frequbigKMha[,6:11]

testKMha<-findOverlaps(frequbigKMha_GRange,lyrGenes)
test2KMha<-pintersect(frequbigKMha_GRange[queryHits(testKMha)],lyrGenes[subjectHits(testKMha)])
test3KMha<-as.data.frame(test2KMha)
frequbigKMha_lyrGenes=mergeByOverlaps(frequbigKMha_GRange,lyrGenes,type=c("any"))
frequbigKMha_lyrGenes_df=data.frame(as.data.frame(frequbigKMha_lyrGenes$frequbigKMha_GRange),as.data.frame(frequbigKMha_lyrGenes$lyrGenes),test3KMha$width)
colnames(frequbigKMha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKMha_lyrGenes_df<-frequbigKMha_lyrGenes_df[(frequbigKMha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKMha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKMha_lyrGenes_df,"KMha_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigKMha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKMha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKMha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup=merge(frequbigKMha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KMha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KMha.xlsx",overwrite=T)

Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKMha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriKowaMiashallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

MKha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MKha.xlsx",1)
KMha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KMha.xlsx",1)

MKha_unique<-MKha[!duplicated(MKha$Lyr_Gene,MKha$Copy_start,MKha$Copy_end),]
KMha_unique<-KMha[!duplicated(KMha$Lyr_Gene,KMha$Copy_start,KMha$Copy_end),]

MKha_unique$Lyr_Gene<-as.character(MKha_unique$Lyr_Gene)
KMha_unique$Lyr_Gene<-as.character(KMha_unique$Lyr_Gene)
MKha_unique$CN_class<-as.character(MKha_unique$CN_class)
KMha_unique$CN_class<-as.character(KMha_unique$CN_class)

MKhaonly<-MKha_unique[!((MKha_unique$Lyr_Gene%in%KMha_unique$Lyr_Gene)&(MKha_unique$CN_class==KMha_unique$CN_class[match(MKha_unique$Lyr_Gene,KMha_unique$Lyr_Gene)])),]
KMhaonly<-KMha_unique[!((KMha_unique$Lyr_Gene%in%MKha_unique$Lyr_Gene)&(KMha_unique$CN_class==MKha_unique$CN_class[match(KMha_unique$Lyr_Gene,MKha_unique$Lyr_Gene)])),]

write.table(MKhaonly,"MKha_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(KMhaonly,"KMha_cnvs_only_strict.table",row.names=F,sep="\t")











