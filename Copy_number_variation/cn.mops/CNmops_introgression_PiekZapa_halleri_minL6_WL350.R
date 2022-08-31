#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesPiekhalleri <- scan("Piek_halleri_full.list",character(),quote="")
BAMFilesZapahalleri <- scan("Zapa_halleri_full.list",character(),quote="")

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

Piekhalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesPiekhalleri,sampleNames=c("Piek_h01","Piek_h03","Piek_h04","Piek_h05comb","Piek_h06comb","Piek_h07","Piek_h10","Piek_h11","Piek_h12","Piek_h13","Piek_h14_2","Piek_h15","Piek_h2"),refSeqName=sort(scaffolds),WL=350) 
Zapahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zako_h01comb","Zako_h03","Zako_h09comb","Zako_h10","Zako_h12","Zapa_h01comb","Zapa_h07","Zapa_h11comb","Zapa_h12"),refSeqName=sort(scaffolds),WL=350) 
PZha <- referencecn.mops(cases=Piekhalleri,controls=Zapahalleri, minWidth = 6, minReadCount=6)
resCNMOPSPZha <- calcIntegerCopyNumbers(PZha)

ZPha <- referencecn.mops(cases=Zapahalleri,controls=Piekhalleri, minWidth = 6, minReadCount=6)
resCNMOPSZPha <- calcIntegerCopyNumbers(ZPha)


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
regionsPZha<-cnvr(resCNMOPSPZha)
singlePZha<-cnvs(resCNMOPSPZha)
testPZha<-findOverlaps(regionsPZha,singlePZha)
test2PZha<- DataFrame(splitAsList(singlePZha$sampleName[subjectHits(testPZha)], queryHits(testPZha)),splitAsList(singlePZha$median[subjectHits(testPZha)], queryHits(testPZha)),splitAsList(singlePZha$mean[subjectHits(testPZha)], queryHits(testPZha)),splitAsList(singlePZha$CN[subjectHits(testPZha)], queryHits(testPZha)))
colnames(test2PZha)<-c("sampleNames","mean","median","CN")
test7PZha<-apply(as.data.frame(test2PZha),1,dupl_cov)
mcols(regionsPZha)<-as.data.frame(as.matrix(test7PZha))

dataPZha<-as.data.frame(regionsPZha)

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

adddataPZha<-do.call("rbind",apply(dataPZha,1,data_change))
dataPZha<-dataPZha[,-6]
frequPZha<-cbind(dataPZha,adddataPZha)

names(frequPZha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3PZha<-frequPZha[count.fields(textConnection(frequPZha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3PZha,"CnvregionswithmeanhalleriPiekZapahalleri_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigPZha<-read.xlsx("CnvregionswithmeanhalleriPiekZapahalleri_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigPZha[,2]<-as.numeric(as.character(frequbigPZha[,2]))
frequbigPZha[,3]<-as.numeric(as.character(frequbigPZha[,3]))

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
frequbigPZha_GRange<-GRanges(seqnames=tolower(frequbigPZha$Scaffold),ranges=IRanges(start=frequbigPZha$start,end=frequbigPZha$end))
values(frequbigPZha_GRange)<-frequbigPZha[,6:11]

testPZha<-findOverlaps(frequbigPZha_GRange,lyrGenes)
test2PZha<-pintersect(frequbigPZha_GRange[queryHits(testPZha)],lyrGenes[subjectHits(testPZha)])
test3PZha<-as.data.frame(test2PZha)
frequbigPZha_lyrGenes=mergeByOverlaps(frequbigPZha_GRange,lyrGenes,type=c("any"))
frequbigPZha_lyrGenes_df=data.frame(as.data.frame(frequbigPZha_lyrGenes$frequbigPZha_GRange),as.data.frame(frequbigPZha_lyrGenes$lyrGenes),test3PZha$width)
colnames(frequbigPZha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigPZha_lyrGenes_df<-frequbigPZha_lyrGenes_df[(frequbigPZha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigPZha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigPZha_lyrGenes_df,"PZha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigPZha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigPZha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigPZha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup=merge(frequbigPZha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PZha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PZha.xlsx",overwrite=T)

Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigPZha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriPiekZapahallericnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Piek halleri#
#################################################

regionsZPha<-cnvr(resCNMOPSZPha)
singleZPha<-cnvs(resCNMOPSZPha)
testZPha<-findOverlaps(regionsZPha,singleZPha)
test2ZPha<- DataFrame(splitAsList(singleZPha$sampleName[subjectHits(testZPha)], queryHits(testZPha)),splitAsList(singleZPha$median[subjectHits(testZPha)], queryHits(testZPha)),splitAsList(singleZPha$mean[subjectHits(testZPha)], queryHits(testZPha)),splitAsList(singleZPha$CN[subjectHits(testZPha)], queryHits(testZPha)))
colnames(test2ZPha)<-c("sampleNames","mean","median","CN")
test7ZPha<-apply(as.data.frame(test2ZPha),1,dupl_cov)
mcols(regionsZPha)<-as.data.frame(as.matrix(test7ZPha))

dataZPha<-as.data.frame(regionsZPha)

adddataZPha<-do.call("rbind",apply(dataZPha,1,data_change))
dataZPha<-dataZPha[,-6]
frequZPha<-cbind(dataZPha,adddataZPha)

names(frequZPha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZPha<-frequZPha[count.fields(textConnection(frequZPha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZPha,"CnvregionswithmeanhalleriZapaPiekhalleri_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZPha<-read.xlsx("CnvregionswithmeanhalleriZapaPiekhalleri_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigZPha[,2]<-as.numeric(as.character(frequbigZPha[,2]))
frequbigZPha[,3]<-as.numeric(as.character(frequbigZPha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZPha_GRange<-GRanges(seqnames=tolower(frequbigZPha$Scaffold),ranges=IRanges(start=frequbigZPha$start,end=frequbigZPha$end))
values(frequbigZPha_GRange)<-frequbigZPha[,6:11]

testZPha<-findOverlaps(frequbigZPha_GRange,lyrGenes)
test2ZPha<-pintersect(frequbigZPha_GRange[queryHits(testZPha)],lyrGenes[subjectHits(testZPha)])
test3ZPha<-as.data.frame(test2ZPha)
frequbigZPha_lyrGenes=mergeByOverlaps(frequbigZPha_GRange,lyrGenes,type=c("any"))
frequbigZPha_lyrGenes_df=data.frame(as.data.frame(frequbigZPha_lyrGenes$frequbigZPha_GRange),as.data.frame(frequbigZPha_lyrGenes$lyrGenes),test3ZPha$width)
colnames(frequbigZPha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZPha_lyrGenes_df<-frequbigZPha_lyrGenes_df[(frequbigZPha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZPha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZPha_lyrGenes_df,"ZPha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZPha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZPha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZPha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup=merge(frequbigZPha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZPha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZPha.xlsx",overwrite=T)

Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZPha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriZapaPiekhallericnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

PZha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PZha.xlsx",1)
ZPha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZPha.xlsx",1)

PZha_unique<-PZha[!duplicated(PZha$Lyr_Gene,PZha$Copy_start,PZha$Copy_end),]
ZPha_unique<-ZPha[!duplicated(ZPha$Lyr_Gene,ZPha$Copy_start,ZPha$Copy_end),]

PZha_unique$Lyr_Gene<-as.character(PZha_unique$Lyr_Gene)
ZPha_unique$Lyr_Gene<-as.character(ZPha_unique$Lyr_Gene)
PZha_unique$CN_class<-as.character(PZha_unique$CN_class)
ZPha_unique$CN_class<-as.character(ZPha_unique$CN_class)

PZhaonly<-PZha_unique[!((PZha_unique$Lyr_Gene%in%ZPha_unique$Lyr_Gene)&(PZha_unique$CN_class==ZPha_unique$CN_class[match(PZha_unique$Lyr_Gene,ZPha_unique$Lyr_Gene)])),]
ZPhaonly<-ZPha_unique[!((ZPha_unique$Lyr_Gene%in%PZha_unique$Lyr_Gene)&(ZPha_unique$CN_class==PZha_unique$CN_class[match(ZPha_unique$Lyr_Gene,PZha_unique$Lyr_Gene)])),]

write.table(PZhaonly,"PZha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(ZPhaonly,"ZPha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











