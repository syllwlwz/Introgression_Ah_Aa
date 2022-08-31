#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesLyrata <- scan("Lyrata_Vinod.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesLyrata)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Lyrata <- getReadCountsFromBAM(BAMFiles=BAMFilesLyrata,sampleNames=c("ICE1","Plech"),refSeqName=sort(scaffolds),WL=350) 
Zapahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zako_h01comb","Zako_h03","Zako_h09comb","Zako_h10","Zako_h12","Zapa_h01comb","Zapa_h07","Zapa_h11comb","Zapa_h12"),refSeqName=sort(scaffolds),WL=350) 
LZha <- referencecn.mops(cases=Lyrata,controls=Zapahalleri, minWidth = 6, minReadCount=6)
resCNMOPSLZha <- calcIntegerCopyNumbers(LZha)

#Lyrataarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesLyrataarenosa,sampleNames=c("Lyrata_h04","Lyrata_h12","Lyrata_h1_7","Lyrata_h21","Lyrata_h22","Lyrata_h23","Lyrata_h24","Lyrata_h25","Lyrata_h26","Lyrata_h27","Lyrata_h28","Lyrata_h29","Lyrata_h30","Lyrata_h31","Lyrata_h32","Lyrata_h33","Lyrata_h35","Lyrata_h36"),refSeqName=sort(scaffolds),WL=500) 
#Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zapa002h13","Zapa002h14","Zapa002h17","Zapa002h18","Zapa002h21","Zapa_h02comb","Zapa_h04comb","Zapa_h07","Zapa_h41","Zapa_h42","Zapa_h43"),refSeqName=sort(scaffolds),WL=500) 
#BKa <- referencecn.mops(cases=Lyrataarenosa,controls=Zapaarenosa, minWidth = 2, minReadCount=6)
#resCNMOPSLZha <- calcIntegerCopyNumbers(BKa)


ZLha <- referencecn.mops(cases=Zapahalleri,controls=Lyrata, minWidth = 6, minReadCount=6)
resCNMOPSZLha <- calcIntegerCopyNumbers(ZLha)


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
regionsLZha<-cnvr(resCNMOPSLZha)
singleLZha<-cnvs(resCNMOPSLZha)
testLZha<-findOverlaps(regionsLZha,singleLZha)
test2LZha<- DataFrame(splitAsList(singleLZha$sampleName[subjectHits(testLZha)], queryHits(testLZha)),splitAsList(singleLZha$median[subjectHits(testLZha)], queryHits(testLZha)),splitAsList(singleLZha$mean[subjectHits(testLZha)], queryHits(testLZha)),splitAsList(singleLZha$CN[subjectHits(testLZha)], queryHits(testLZha)))
colnames(test2LZha)<-c("sampleNames","mean","median","CN")
test7LZha<-apply(as.data.frame(test2LZha),1,dupl_cov)
mcols(regionsLZha)<-as.data.frame(as.matrix(test7LZha))

dataLZha<-as.data.frame(regionsLZha)

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

adddataLZha<-do.call("rbind",apply(dataLZha,1,data_change))
dataLZha<-dataLZha[,-6]
frequLZha<-cbind(dataLZha,adddataLZha)

names(frequLZha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3LZha<-frequLZha[count.fields(textConnection(frequLZha$sampleNames),sep=",")>=2,]

require(openxlsx)
write.xlsx(frequ3LZha,"CnvregionswithmeanhalleriLyrataZapahalleri_minRead6_WL350_minL6_min2samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigLZha<-read.xlsx("CnvregionswithmeanhalleriLyrataZapahalleri_minRead6_WL350_minL6_min2samples.xlsx",1)
frequbigLZha[,2]<-as.numeric(as.character(frequbigLZha[,2]))
frequbigLZha[,3]<-as.numeric(as.character(frequbigLZha[,3]))

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
frequbigLZha_GRange<-GRanges(seqnames=tolower(frequbigLZha$Scaffold),ranges=IRanges(start=frequbigLZha$start,end=frequbigLZha$end))
values(frequbigLZha_GRange)<-frequbigLZha[,6:11]

testLZha<-findOverlaps(frequbigLZha_GRange,lyrGenes)
test2LZha<-pintersect(frequbigLZha_GRange[queryHits(testLZha)],lyrGenes[subjectHits(testLZha)])
test3LZha<-as.data.frame(test2LZha)
frequbigLZha_lyrGenes=mergeByOverlaps(frequbigLZha_GRange,lyrGenes,type=c("any"))
frequbigLZha_lyrGenes_df=data.frame(as.data.frame(frequbigLZha_lyrGenes$frequbigLZha_GRange),as.data.frame(frequbigLZha_lyrGenes$lyrGenes),test3LZha$width)
colnames(frequbigLZha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigLZha_lyrGenes_df<-frequbigLZha_lyrGenes_df[(frequbigLZha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigLZha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigLZha_lyrGenes_df,"LZha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigLZha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigLZha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigLZha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup=merge(frequbigLZha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LZha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LZha.xlsx",overwrite=T)

Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigLZha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriLyrataZapahallericnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Lyrata halleri#
#################################################

regionsZLha<-cnvr(resCNMOPSZLha)
singleZLha<-cnvs(resCNMOPSZLha)
testZLha<-findOverlaps(regionsZLha,singleZLha)
test2ZLha<- DataFrame(splitAsList(singleZLha$sampleName[subjectHits(testZLha)], queryHits(testZLha)),splitAsList(singleZLha$median[subjectHits(testZLha)], queryHits(testZLha)),splitAsList(singleZLha$mean[subjectHits(testZLha)], queryHits(testZLha)),splitAsList(singleZLha$CN[subjectHits(testZLha)], queryHits(testZLha)))
colnames(test2ZLha)<-c("sampleNames","mean","median","CN")
test7ZLha<-apply(as.data.frame(test2ZLha),1,dupl_cov)
mcols(regionsZLha)<-as.data.frame(as.matrix(test7ZLha))

dataZLha<-as.data.frame(regionsZLha)

adddataZLha<-do.call("rbind",apply(dataZLha,1,data_change))
dataZLha<-dataZLha[,-6]
frequZLha<-cbind(dataZLha,adddataZLha)

names(frequZLha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZLha<-frequZLha[count.fields(textConnection(frequZLha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZLha,"CnvregionswithmeanhalleriZapaLyrata_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZLha<-read.xlsx("CnvregionswithmeanhalleriZapaLyrata_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigZLha[,2]<-as.numeric(as.character(frequbigZLha[,2]))
frequbigZLha[,3]<-as.numeric(as.character(frequbigZLha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZLha_GRange<-GRanges(seqnames=tolower(frequbigZLha$Scaffold),ranges=IRanges(start=frequbigZLha$start,end=frequbigZLha$end))
values(frequbigZLha_GRange)<-frequbigZLha[,6:11]

testZLha<-findOverlaps(frequbigZLha_GRange,lyrGenes)
test2ZLha<-pintersect(frequbigZLha_GRange[queryHits(testZLha)],lyrGenes[subjectHits(testZLha)])
test3ZLha<-as.data.frame(test2ZLha)
frequbigZLha_lyrGenes=mergeByOverlaps(frequbigZLha_GRange,lyrGenes,type=c("any"))
frequbigZLha_lyrGenes_df=data.frame(as.data.frame(frequbigZLha_lyrGenes$frequbigZLha_GRange),as.data.frame(frequbigZLha_lyrGenes$lyrGenes),test3ZLha$width)
colnames(frequbigZLha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZLha_lyrGenes_df<-frequbigZLha_lyrGenes_df[(frequbigZLha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZLha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZLha_lyrGenes_df,"ZLha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZLha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZLha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZLha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup=merge(frequbigZLha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZLha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZLha.xlsx",overwrite=T)

Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZLha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriZapaLyratacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

LZha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LZha.xlsx",1)
ZLha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZLha.xlsx",1)

LZha_unique<-LZha[!duplicated(LZha$Lyr_Gene,LZha$Copy_start,LZha$Copy_end),]
ZLha_unique<-ZLha[!duplicated(ZLha$Lyr_Gene,ZLha$Copy_start,ZLha$Copy_end),]

LZha_unique$Lyr_Gene<-as.character(LZha_unique$Lyr_Gene)
ZLha_unique$Lyr_Gene<-as.character(ZLha_unique$Lyr_Gene)
LZha_unique$CN_class<-as.character(LZha_unique$CN_class)
ZLha_unique$CN_class<-as.character(ZLha_unique$CN_class)

LZhaonly<-LZha_unique[!((LZha_unique$Lyr_Gene%in%ZLha_unique$Lyr_Gene)&(LZha_unique$CN_class==ZLha_unique$CN_class[match(LZha_unique$Lyr_Gene,ZLha_unique$Lyr_Gene)])),]
ZLhaonly<-ZLha_unique[!((ZLha_unique$Lyr_Gene%in%LZha_unique$Lyr_Gene)&(ZLha_unique$CN_class==LZha_unique$CN_class[match(ZLha_unique$Lyr_Gene,LZha_unique$Lyr_Gene)])),]

write.table(LZhaonly,"LZha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(ZLhaonly,"ZLha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











