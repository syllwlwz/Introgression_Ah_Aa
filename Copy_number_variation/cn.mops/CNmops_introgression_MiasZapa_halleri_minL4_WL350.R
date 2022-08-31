#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesMiashalleri <- scan("Mias_halleri_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesMiashalleri)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Miashalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesMiashalleri,sampleNames=c("Mias002h19","Mias002h20","Mias004h02","Mias004h06","Mias004h07","Mias004h08","Mias009h03","Mias009h05","Mias_a54","Mias_h09","Mias_h41","Mias_h43","Mias_h44","Mias_h45","Mias_h46"),refSeqName=sort(scaffolds),WL=350) 
Zapahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zako_h01comb","Zako_h03","Zako_h09comb","Zako_h10","Zako_h12","Zapa_h01comb","Zapa_h07","Zapa_h11comb","Zapa_h12"),refSeqName=sort(scaffolds),WL=350) 
MZha <- referencecn.mops(cases=Miashalleri,controls=Zapahalleri, minWidth = 4, minReadCount=6)
resCNMOPSMZha <- calcIntegerCopyNumbers(MZha)

ZMha <- referencecn.mops(cases=Zapahalleri,controls=Miashalleri, minWidth = 4, minReadCount=6)
resCNMOPSZMha <- calcIntegerCopyNumbers(ZMha)


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
regionsMZha<-cnvr(resCNMOPSMZha)
singleMZha<-cnvs(resCNMOPSMZha)
testMZha<-findOverlaps(regionsMZha,singleMZha)
test2MZha<- DataFrame(splitAsList(singleMZha$sampleName[subjectHits(testMZha)], queryHits(testMZha)),splitAsList(singleMZha$median[subjectHits(testMZha)], queryHits(testMZha)),splitAsList(singleMZha$mean[subjectHits(testMZha)], queryHits(testMZha)),splitAsList(singleMZha$CN[subjectHits(testMZha)], queryHits(testMZha)))
colnames(test2MZha)<-c("sampleNames","mean","median","CN")
test7MZha<-apply(as.data.frame(test2MZha),1,dupl_cov)
mcols(regionsMZha)<-as.data.frame(as.matrix(test7MZha))

dataMZha<-as.data.frame(regionsMZha)

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

adddataMZha<-do.call("rbind",apply(dataMZha,1,data_change))
dataMZha<-dataMZha[,-6]
frequMZha<-cbind(dataMZha,adddataMZha)

names(frequMZha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3MZha<-frequMZha[count.fields(textConnection(frequMZha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3MZha,"CnvregionswithmeanhalleriMiasZapahalleri_minRead6_WL350_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigMZha<-read.xlsx("CnvregionswithmeanhalleriMiasZapahalleri_minRead6_WL350_minL4_min5samples.xlsx",1)
frequbigMZha[,2]<-as.numeric(as.character(frequbigMZha[,2]))
frequbigMZha[,3]<-as.numeric(as.character(frequbigMZha[,3]))

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
frequbigMZha_GRange<-GRanges(seqnames=tolower(frequbigMZha$Scaffold),ranges=IRanges(start=frequbigMZha$start,end=frequbigMZha$end))
values(frequbigMZha_GRange)<-frequbigMZha[,6:11]

testMZha<-findOverlaps(frequbigMZha_GRange,lyrGenes)
test2MZha<-pintersect(frequbigMZha_GRange[queryHits(testMZha)],lyrGenes[subjectHits(testMZha)])
test3MZha<-as.data.frame(test2MZha)
frequbigMZha_lyrGenes=mergeByOverlaps(frequbigMZha_GRange,lyrGenes,type=c("any"))
frequbigMZha_lyrGenes_df=data.frame(as.data.frame(frequbigMZha_lyrGenes$frequbigMZha_GRange),as.data.frame(frequbigMZha_lyrGenes$lyrGenes),test3MZha$width)
colnames(frequbigMZha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigMZha_lyrGenes_df<-frequbigMZha_lyrGenes_df[(frequbigMZha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigMZha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigMZha_lyrGenes_df,"MZha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigMZha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigMZha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigMZha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup=merge(frequbigMZha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MZha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MZha.xlsx",overwrite=T)

Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigMZha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriMiasZapahallericnvregions_minRead6_WL350_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Mias halleri#
#################################################

regionsZMha<-cnvr(resCNMOPSZMha)
singleZMha<-cnvs(resCNMOPSZMha)
testZMha<-findOverlaps(regionsZMha,singleZMha)
test2ZMha<- DataFrame(splitAsList(singleZMha$sampleName[subjectHits(testZMha)], queryHits(testZMha)),splitAsList(singleZMha$median[subjectHits(testZMha)], queryHits(testZMha)),splitAsList(singleZMha$mean[subjectHits(testZMha)], queryHits(testZMha)),splitAsList(singleZMha$CN[subjectHits(testZMha)], queryHits(testZMha)))
colnames(test2ZMha)<-c("sampleNames","mean","median","CN")
test7ZMha<-apply(as.data.frame(test2ZMha),1,dupl_cov)
mcols(regionsZMha)<-as.data.frame(as.matrix(test7ZMha))

dataZMha<-as.data.frame(regionsZMha)

adddataZMha<-do.call("rbind",apply(dataZMha,1,data_change))
dataZMha<-dataZMha[,-6]
frequZMha<-cbind(dataZMha,adddataZMha)

names(frequZMha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZMha<-frequZMha[count.fields(textConnection(frequZMha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZMha,"CnvregionswithmeanhalleriZapaMiashalleri_minRead6_WL350_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZMha<-read.xlsx("CnvregionswithmeanhalleriZapaMiashalleri_minRead6_WL350_minL4_min5samples.xlsx",1)
frequbigZMha[,2]<-as.numeric(as.character(frequbigZMha[,2]))
frequbigZMha[,3]<-as.numeric(as.character(frequbigZMha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZMha_GRange<-GRanges(seqnames=tolower(frequbigZMha$Scaffold),ranges=IRanges(start=frequbigZMha$start,end=frequbigZMha$end))
values(frequbigZMha_GRange)<-frequbigZMha[,6:11]

testZMha<-findOverlaps(frequbigZMha_GRange,lyrGenes)
test2ZMha<-pintersect(frequbigZMha_GRange[queryHits(testZMha)],lyrGenes[subjectHits(testZMha)])
test3ZMha<-as.data.frame(test2ZMha)
frequbigZMha_lyrGenes=mergeByOverlaps(frequbigZMha_GRange,lyrGenes,type=c("any"))
frequbigZMha_lyrGenes_df=data.frame(as.data.frame(frequbigZMha_lyrGenes$frequbigZMha_GRange),as.data.frame(frequbigZMha_lyrGenes$lyrGenes),test3ZMha$width)
colnames(frequbigZMha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZMha_lyrGenes_df<-frequbigZMha_lyrGenes_df[(frequbigZMha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZMha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZMha_lyrGenes_df,"ZMha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZMha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZMha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZMha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup=merge(frequbigZMha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZMha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZMha.xlsx",overwrite=T)

Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZMha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriZapaMiashallericnvregions_minRead6_WL350_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

MZha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MZha.xlsx",1)
ZMha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZMha.xlsx",1)

MZha_unique<-MZha[!duplicated(MZha$Lyr_Gene,MZha$Copy_start,MZha$Copy_end),]
ZMha_unique<-ZMha[!duplicated(ZMha$Lyr_Gene,ZMha$Copy_start,ZMha$Copy_end),]

MZha_unique$Lyr_Gene<-as.character(MZha_unique$Lyr_Gene)
ZMha_unique$Lyr_Gene<-as.character(ZMha_unique$Lyr_Gene)
MZha_unique$CN_class<-as.character(MZha_unique$CN_class)
ZMha_unique$CN_class<-as.character(ZMha_unique$CN_class)

MZhaonly<-MZha_unique[!((MZha_unique$Lyr_Gene%in%ZMha_unique$Lyr_Gene)&(MZha_unique$CN_class==ZMha_unique$CN_class[match(MZha_unique$Lyr_Gene,ZMha_unique$Lyr_Gene)])),]
ZMhaonly<-ZMha_unique[!((ZMha_unique$Lyr_Gene%in%MZha_unique$Lyr_Gene)&(ZMha_unique$CN_class==MZha_unique$CN_class[match(ZMha_unique$Lyr_Gene,MZha_unique$Lyr_Gene)])),]

write.table(MZhaonly,"MZha_cnvs_only_strict_minL4_WL350.table",row.names=F,sep="\t")
write.table(ZMhaonly,"ZMha_cnvs_only_strict_minL4_WL350.table",row.names=F,sep="\t")











