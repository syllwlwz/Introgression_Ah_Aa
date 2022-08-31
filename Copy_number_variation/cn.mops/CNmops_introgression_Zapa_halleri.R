#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesZapahalleri <- scan("Zapa_halleri_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesZapahalleri)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Zapahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zako_h01comb","Zako_h03","Zako_h09comb","Zako_h10","Zako_h12","Zapa_h01comb","Zapa_h07","Zapa_h11comb","Zapa_h12"),refSeqName=sort(scaffolds),WL=500) 
Kowahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesKowahalleri,sampleNames=c("Kowa002h13","Kowa002h14","Kowa002h17","Kowa002h18","Kowa002h21","Kowa_h02comb","Kowa_h04comb","Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43"),refSeqName=sort(scaffolds),WL=500) 
ZKha <- referencecn.mops(cases=Zapahalleri,controls=Kowahalleri, minWidth = 2, minReadCount=6)
resCNMOPSZKha <- calcIntegerCopyNumbers(ZKha)

KZha <- referencecn.mops(cases=Kowahalleri,controls=Zapahalleri, minWidth = 2, minReadCount=6)
resCNMOPSKZha <- calcIntegerCopyNumbers(KZha)


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
regionsZKha<-cnvr(resCNMOPSZKha)
singleZKha<-cnvs(resCNMOPSZKha)
testZKha<-findOverlaps(regionsZKha,singleZKha)
test2ZKha<- DataFrame(splitAsList(singleZKha$sampleName[subjectHits(testZKha)], queryHits(testZKha)),splitAsList(singleZKha$median[subjectHits(testZKha)], queryHits(testZKha)),splitAsList(singleZKha$mean[subjectHits(testZKha)], queryHits(testZKha)),splitAsList(singleZKha$CN[subjectHits(testZKha)], queryHits(testZKha)))
colnames(test2ZKha)<-c("sampleNames","mean","median","CN")
test7ZKha<-apply(as.data.frame(test2ZKha),1,dupl_cov)
mcols(regionsZKha)<-as.data.frame(as.matrix(test7ZKha))

dataZKha<-as.data.frame(regionsZKha)

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

adddataZKha<-do.call("rbind",apply(dataZKha,1,data_change))
dataZKha<-dataZKha[,-6]
frequZKha<-cbind(dataZKha,adddataZKha)

names(frequZKha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZKha<-frequZKha[count.fields(textConnection(frequZKha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZKha,"CnvregionswithmeanhalleriZapaKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZKha<-read.xlsx("CnvregionswithmeanhalleriZapaKowahalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigZKha[,2]<-as.numeric(as.character(frequbigZKha[,2]))
frequbigZKha[,3]<-as.numeric(as.character(frequbigZKha[,3]))

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
frequbigZKha_GRange<-GRanges(seqnames=tolower(frequbigZKha$Scaffold),ranges=IRanges(start=frequbigZKha$start,end=frequbigZKha$end))
values(frequbigZKha_GRange)<-frequbigZKha[,6:11]

testZKha<-findOverlaps(frequbigZKha_GRange,lyrGenes)
test2ZKha<-pintersect(frequbigZKha_GRange[queryHits(testZKha)],lyrGenes[subjectHits(testZKha)])
test3ZKha<-as.data.frame(test2ZKha)
frequbigZKha_lyrGenes=mergeByOverlaps(frequbigZKha_GRange,lyrGenes,type=c("any"))
frequbigZKha_lyrGenes_df=data.frame(as.data.frame(frequbigZKha_lyrGenes$frequbigZKha_GRange),as.data.frame(frequbigZKha_lyrGenes$lyrGenes),test3ZKha$width)
colnames(frequbigZKha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZKha_lyrGenes_df<-frequbigZKha_lyrGenes_df[(frequbigZKha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZKha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZKha_lyrGenes_df,"ZKoha_temp.xlsx",col.names=TRUE,row.names=FALSE)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigZKha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZKha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZKha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup=merge(frequbigZKha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZKoha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZKoha.xlsx",overwrite=T)

Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigZKha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriZapaKowahallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Zapa halleri#
#################################################

regionsKZha<-cnvr(resCNMOPSKZha)
singleKZha<-cnvs(resCNMOPSKZha)
testKZha<-findOverlaps(regionsKZha,singleKZha)
test2KZha<- DataFrame(splitAsList(singleKZha$sampleName[subjectHits(testKZha)], queryHits(testKZha)),splitAsList(singleKZha$median[subjectHits(testKZha)], queryHits(testKZha)),splitAsList(singleKZha$mean[subjectHits(testKZha)], queryHits(testKZha)),splitAsList(singleKZha$CN[subjectHits(testKZha)], queryHits(testKZha)))
colnames(test2KZha)<-c("sampleNames","mean","median","CN")
test7KZha<-apply(as.data.frame(test2KZha),1,dupl_cov)
mcols(regionsKZha)<-as.data.frame(as.matrix(test7KZha))

dataKZha<-as.data.frame(regionsKZha)

adddataKZha<-do.call("rbind",apply(dataKZha,1,data_change))
dataKZha<-dataKZha[,-6]
frequKZha<-cbind(dataKZha,adddataKZha)

names(frequKZha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KZha<-frequKZha[count.fields(textConnection(frequKZha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KZha,"CnvregionswithmeanhalleriKowaZapahalleri_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKZha<-read.xlsx("CnvregionswithmeanhalleriKowaZapahalleri_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigKZha[,2]<-as.numeric(as.character(frequbigKZha[,2]))
frequbigKZha[,3]<-as.numeric(as.character(frequbigKZha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKZha_GRange<-GRanges(seqnames=tolower(frequbigKZha$Scaffold),ranges=IRanges(start=frequbigKZha$start,end=frequbigKZha$end))
values(frequbigKZha_GRange)<-frequbigKZha[,6:11]

testKZha<-findOverlaps(frequbigKZha_GRange,lyrGenes)
test2KZha<-pintersect(frequbigKZha_GRange[queryHits(testKZha)],lyrGenes[subjectHits(testKZha)])
test3KZha<-as.data.frame(test2KZha)
frequbigKZha_lyrGenes=mergeByOverlaps(frequbigKZha_GRange,lyrGenes,type=c("any"))
frequbigKZha_lyrGenes_df=data.frame(as.data.frame(frequbigKZha_lyrGenes$frequbigKZha_GRange),as.data.frame(frequbigKZha_lyrGenes$lyrGenes),test3KZha$width)
colnames(frequbigKZha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKZha_lyrGenes_df<-frequbigKZha_lyrGenes_df[(frequbigKZha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKZha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKZha_lyrGenes_df,"KoZha_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigKZha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKZha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKZha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup=merge(frequbigKZha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KoZha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KoZha.xlsx",overwrite=T)

Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKZha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriKowaZapahallericnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

ZKha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZKoha.xlsx",1)
KZha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KoZha.xlsx",1)

ZKha_unique<-ZKha[!duplicated(ZKha$Lyr_Gene,ZKha$Copy_start,ZKha$Copy_end),]
KZha_unique<-KZha[!duplicated(KZha$Lyr_Gene,KZha$Copy_start,KZha$Copy_end),]

ZKha_unique$Lyr_Gene<-as.character(ZKha_unique$Lyr_Gene)
KZha_unique$Lyr_Gene<-as.character(KZha_unique$Lyr_Gene)
ZKha_unique$CN_class<-as.character(ZKha_unique$CN_class)
KZha_unique$CN_class<-as.character(KZha_unique$CN_class)

ZKhaonly<-ZKha_unique[!((ZKha_unique$Lyr_Gene%in%KZha_unique$Lyr_Gene)&(ZKha_unique$CN_class==KZha_unique$CN_class[match(ZKha_unique$Lyr_Gene,KZha_unique$Lyr_Gene)])),]
KZhaonly<-KZha_unique[!((KZha_unique$Lyr_Gene%in%ZKha_unique$Lyr_Gene)&(KZha_unique$CN_class==ZKha_unique$CN_class[match(KZha_unique$Lyr_Gene,ZKha_unique$Lyr_Gene)])),]

write.table(ZKhaonly,"ZKoha_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(KZhaonly,"KoZha_cnvs_only_strict.table",row.names=F,sep="\t")











