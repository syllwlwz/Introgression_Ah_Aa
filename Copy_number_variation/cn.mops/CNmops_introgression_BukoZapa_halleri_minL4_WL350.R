#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesBukohalleri <- scan("Buko_halleri_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesBukohalleri)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Bukohalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesBukohalleri,sampleNames=c("Buko_h04","Buko_h12","Buko_h1_7","Buko_h21","Buko_h22","Buko_h23","Buko_h24","Buko_h25","Buko_h26","Buko_h27","Buko_h28","Buko_h29","Buko_h30","Buko_h31","Buko_h32","Buko_h33","Buko_h35","Buko_h36"),refSeqName=sort(scaffolds),WL=350) 
Zapahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zako_h01comb","Zako_h03","Zako_h09comb","Zako_h10","Zako_h12","Zapa_h01comb","Zapa_h07","Zapa_h11comb","Zapa_h12"),refSeqName=sort(scaffolds),WL=350) 
BZha <- referencecn.mops(cases=Bukohalleri,controls=Zapahalleri, minWidth = 4, minReadCount=6)
resCNMOPSBZha <- calcIntegerCopyNumbers(BZha)

#Bukoarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesBukoarenosa,sampleNames=c("Buko_h04","Buko_h12","Buko_h1_7","Buko_h21","Buko_h22","Buko_h23","Buko_h24","Buko_h25","Buko_h26","Buko_h27","Buko_h28","Buko_h29","Buko_h30","Buko_h31","Buko_h32","Buko_h33","Buko_h35","Buko_h36"),refSeqName=sort(scaffolds),WL=500) 
#Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zapa002h13","Zapa002h14","Zapa002h17","Zapa002h18","Zapa002h21","Zapa_h02comb","Zapa_h04comb","Zapa_h07","Zapa_h41","Zapa_h42","Zapa_h43"),refSeqName=sort(scaffolds),WL=500) 
#BKa <- referencecn.mops(cases=Bukoarenosa,controls=Zapaarenosa, minWidth = 2, minReadCount=6)
#resCNMOPSBZha <- calcIntegerCopyNumbers(BKa)


ZBha <- referencecn.mops(cases=Zapahalleri,controls=Bukohalleri, minWidth = 4, minReadCount=6)
resCNMOPSZBha <- calcIntegerCopyNumbers(ZBha)


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
regionsBZha<-cnvr(resCNMOPSBZha)
singleBZha<-cnvs(resCNMOPSBZha)
testBZha<-findOverlaps(regionsBZha,singleBZha)
test2BZha<- DataFrame(splitAsList(singleBZha$sampleName[subjectHits(testBZha)], queryHits(testBZha)),splitAsList(singleBZha$median[subjectHits(testBZha)], queryHits(testBZha)),splitAsList(singleBZha$mean[subjectHits(testBZha)], queryHits(testBZha)),splitAsList(singleBZha$CN[subjectHits(testBZha)], queryHits(testBZha)))
colnames(test2BZha)<-c("sampleNames","mean","median","CN")
test7BZha<-apply(as.data.frame(test2BZha),1,dupl_cov)
mcols(regionsBZha)<-as.data.frame(as.matrix(test7BZha))

dataBZha<-as.data.frame(regionsBZha)

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

adddataBZha<-do.call("rbind",apply(dataBZha,1,data_change))
dataBZha<-dataBZha[,-6]
frequBZha<-cbind(dataBZha,adddataBZha)

names(frequBZha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3BZha<-frequBZha[count.fields(textConnection(frequBZha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3BZha,"CnvregionswithmeanhalleriBukoZapahalleri_minRead6_WL350_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigBZha<-read.xlsx("CnvregionswithmeanhalleriBukoZapahalleri_minRead6_WL350_minL4_min5samples.xlsx",1)
frequbigBZha[,2]<-as.numeric(as.character(frequbigBZha[,2]))
frequbigBZha[,3]<-as.numeric(as.character(frequbigBZha[,3]))

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
frequbigBZha_GRange<-GRanges(seqnames=tolower(frequbigBZha$Scaffold),ranges=IRanges(start=frequbigBZha$start,end=frequbigBZha$end))
values(frequbigBZha_GRange)<-frequbigBZha[,6:11]

testBZha<-findOverlaps(frequbigBZha_GRange,lyrGenes)
test2BZha<-pintersect(frequbigBZha_GRange[queryHits(testBZha)],lyrGenes[subjectHits(testBZha)])
test3BZha<-as.data.frame(test2BZha)
frequbigBZha_lyrGenes=mergeByOverlaps(frequbigBZha_GRange,lyrGenes,type=c("any"))
frequbigBZha_lyrGenes_df=data.frame(as.data.frame(frequbigBZha_lyrGenes$frequbigBZha_GRange),as.data.frame(frequbigBZha_lyrGenes$lyrGenes),test3BZha$width)
colnames(frequbigBZha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigBZha_lyrGenes_df<-frequbigBZha_lyrGenes_df[(frequbigBZha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigBZha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigBZha_lyrGenes_df,"BZha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigBZha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigBZha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigBZha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup=merge(frequbigBZha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BZha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BZha.xlsx",overwrite=T)

Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigBZha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriBukoZapahallericnvregions_minRead6_WL350_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Buko halleri#
#################################################

regionsZBha<-cnvr(resCNMOPSZBha)
singleZBha<-cnvs(resCNMOPSZBha)
testZBha<-findOverlaps(regionsZBha,singleZBha)
test2ZBha<- DataFrame(splitAsList(singleZBha$sampleName[subjectHits(testZBha)], queryHits(testZBha)),splitAsList(singleZBha$median[subjectHits(testZBha)], queryHits(testZBha)),splitAsList(singleZBha$mean[subjectHits(testZBha)], queryHits(testZBha)),splitAsList(singleZBha$CN[subjectHits(testZBha)], queryHits(testZBha)))
colnames(test2ZBha)<-c("sampleNames","mean","median","CN")
test7ZBha<-apply(as.data.frame(test2ZBha),1,dupl_cov)
mcols(regionsZBha)<-as.data.frame(as.matrix(test7ZBha))

dataZBha<-as.data.frame(regionsZBha)

adddataZBha<-do.call("rbind",apply(dataZBha,1,data_change))
dataZBha<-dataZBha[,-6]
frequZBha<-cbind(dataZBha,adddataZBha)

names(frequZBha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZBha<-frequZBha[count.fields(textConnection(frequZBha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZBha,"CnvregionswithmeanhalleriZapaBukohalleri_minRead6_WL350_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZBha<-read.xlsx("CnvregionswithmeanhalleriZapaBukohalleri_minRead6_WL350_minL4_min5samples.xlsx",1)
frequbigZBha[,2]<-as.numeric(as.character(frequbigZBha[,2]))
frequbigZBha[,3]<-as.numeric(as.character(frequbigZBha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZBha_GRange<-GRanges(seqnames=tolower(frequbigZBha$Scaffold),ranges=IRanges(start=frequbigZBha$start,end=frequbigZBha$end))
values(frequbigZBha_GRange)<-frequbigZBha[,6:11]

testZBha<-findOverlaps(frequbigZBha_GRange,lyrGenes)
test2ZBha<-pintersect(frequbigZBha_GRange[queryHits(testZBha)],lyrGenes[subjectHits(testZBha)])
test3ZBha<-as.data.frame(test2ZBha)
frequbigZBha_lyrGenes=mergeByOverlaps(frequbigZBha_GRange,lyrGenes,type=c("any"))
frequbigZBha_lyrGenes_df=data.frame(as.data.frame(frequbigZBha_lyrGenes$frequbigZBha_GRange),as.data.frame(frequbigZBha_lyrGenes$lyrGenes),test3ZBha$width)
colnames(frequbigZBha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZBha_lyrGenes_df<-frequbigZBha_lyrGenes_df[(frequbigZBha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZBha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZBha_lyrGenes_df,"ZBha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZBha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZBha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZBha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup=merge(frequbigZBha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZBha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZBha.xlsx",overwrite=T)

Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZBha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriZapaBukohallericnvregions_minRead6_WL350_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GeneshalleriWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

BZha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BZha.xlsx",1)
ZBha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZBha.xlsx",1)

BZha_unique<-BZha[!duplicated(BZha$Lyr_Gene,BZha$Copy_start,BZha$Copy_end),]
ZBha_unique<-ZBha[!duplicated(ZBha$Lyr_Gene,ZBha$Copy_start,ZBha$Copy_end),]

BZha_unique$Lyr_Gene<-as.character(BZha_unique$Lyr_Gene)
ZBha_unique$Lyr_Gene<-as.character(ZBha_unique$Lyr_Gene)
BZha_unique$CN_class<-as.character(BZha_unique$CN_class)
ZBha_unique$CN_class<-as.character(ZBha_unique$CN_class)

BZhaonly<-BZha_unique[!((BZha_unique$Lyr_Gene%in%ZBha_unique$Lyr_Gene)&(BZha_unique$CN_class==ZBha_unique$CN_class[match(BZha_unique$Lyr_Gene,ZBha_unique$Lyr_Gene)])),]
ZBhaonly<-ZBha_unique[!((ZBha_unique$Lyr_Gene%in%BZha_unique$Lyr_Gene)&(ZBha_unique$CN_class==BZha_unique$CN_class[match(ZBha_unique$Lyr_Gene,BZha_unique$Lyr_Gene)])),]

write.table(BZhaonly,"BZha_cnvs_only_strict_minL4_WL350.table",row.names=F,sep="\t")
write.table(ZBhaonly,"ZBha_cnvs_only_strict_minL4_WL350.table",row.names=F,sep="\t")











