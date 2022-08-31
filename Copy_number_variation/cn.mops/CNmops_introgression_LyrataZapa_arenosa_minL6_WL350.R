#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesLyrata <- scan("Lyrata_Vinod.list",character(),quote="")
BAMFilesZapaarenosa <- scan("Zapa_arenosa_full.list",character(),quote="")

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
Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapaarenosa,sampleNames=c("Zapa002a01","Zapa002a02","Zapa002a03","Zapa002a04","Zapa002a05","Zapa002a07","Zapa004a11","Zapa008a19","Zapa_11a03x12a01_l","Zapa_a09xa03_b"),refSeqName=sort(scaffolds),WL=350)
LZar <- referencecn.mops(cases=Lyrata,controls=Zapaarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSLZar <- calcIntegerCopyNumbers(LZar)

ZLyr <- referencecn.mops(cases=Zapaarenosa,controls=Lyrata, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSZLyr <- calcIntegerCopyNumbers(ZLyr)


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
regionsLZar<-cnvr(resCNMOPSLZar)
singleLZar<-cnvs(resCNMOPSLZar)
testLZar<-findOverlaps(regionsLZar,singleLZar)
test2LZar<- DataFrame(splitAsList(singleLZar$sampleName[subjectHits(testLZar)], queryHits(testLZar)),splitAsList(singleLZar$median[subjectHits(testLZar)], queryHits(testLZar)),splitAsList(singleLZar$mean[subjectHits(testLZar)], queryHits(testLZar)),splitAsList(singleLZar$CN[subjectHits(testLZar)], queryHits(testLZar)))
colnames(test2LZar)<-c("sampleNames","mean","median","CN")
test7LZar<-apply(as.data.frame(test2LZar),1,dupl_cov)
mcols(regionsLZar)<-as.data.frame(as.matrix(test7LZar))

dataLZar<-as.data.frame(regionsLZar)

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

adddataLZar<-do.call("rbind",apply(dataLZar,1,data_change))
dataLZar<-dataLZar[,-6]
frequLZar<-cbind(dataLZar,adddataLZar)

names(frequLZar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3LZar<-frequLZar[count.fields(textConnection(frequLZar$sampleNames),sep=",")>=2,]

require(openxlsx)
write.xlsx(frequ3LZar,"CnvregionswithmeanarenosaLyrataZapaarenosa_minRead6_WL350_minL6_min2samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
write.xlsx(frequLZar,"CnvregionswithmeanarenosaLyrataZapaarenosa_minRead6_WL350_minL6_test.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigLZar<-read.xlsx("CnvregionswithmeanarenosaLyrataZapaarenosa_minRead6_WL350_minL6_min2samples.xlsx",1)
frequbigLZar[,2]<-as.numeric(as.character(frequbigLZar[,2]))
frequbigLZar[,3]<-as.numeric(as.character(frequbigLZar[,3]))

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
frequbigLZar_GRange<-GRanges(seqnames=tolower(frequbigLZar$Scaffold),ranges=IRanges(start=frequbigLZar$start,end=frequbigLZar$end))
values(frequbigLZar_GRange)<-frequbigLZar[,6:11]

testLZar<-findOverlaps(frequbigLZar_GRange,lyrGenes)
test2LZar<-pintersect(frequbigLZar_GRange[queryHits(testLZar)],lyrGenes[subjectHits(testLZar)])
test3LZar<-as.data.frame(test2LZar)
frequbigLZar_lyrGenes=mergeByOverlaps(frequbigLZar_GRange,lyrGenes,type=c("any"))
frequbigLZar_lyrGenes_df=data.frame(as.data.frame(frequbigLZar_lyrGenes$frequbigLZar_GRange),as.data.frame(frequbigLZar_lyrGenes$lyrGenes),test3LZar$width)
colnames(frequbigLZar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigLZar_lyrGenes_df<-frequbigLZar_lyrGenes_df[(frequbigLZar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigLZar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigLZar_lyrGenes_df,"LZar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigLZar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigLZar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigLZar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup=merge(frequbigLZar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LZar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LZar.xlsx",overwrite=T)

Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigLZar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaLyrataZapaarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Genesarenosacnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Lyrata arenosa#
#################################################

regionsZLyr<-cnvr(resCNMOPSZLyr)
singleZLyr<-cnvs(resCNMOPSZLyr)
testZLyr<-findOverlaps(regionsZLyr,singleZLyr)
test2ZLyr<- DataFrame(splitAsList(singleZLyr$sampleName[subjectHits(testZLyr)], queryHits(testZLyr)),splitAsList(singleZLyr$median[subjectHits(testZLyr)], queryHits(testZLyr)),splitAsList(singleZLyr$mean[subjectHits(testZLyr)], queryHits(testZLyr)),splitAsList(singleZLyr$CN[subjectHits(testZLyr)], queryHits(testZLyr)))
colnames(test2ZLyr)<-c("sampleNames","mean","median","CN")
test7ZLyr<-apply(as.data.frame(test2ZLyr),1,dupl_cov)
mcols(regionsZLyr)<-as.data.frame(as.matrix(test7ZLyr))

dataZLyr<-as.data.frame(regionsZLyr)

adddataZLyr<-do.call("rbind",apply(dataZLyr,1,data_change))
dataZLyr<-dataZLyr[,-6]
frequZLyr<-cbind(dataZLyr,adddataZLyr)

names(frequZLyr)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZLyr<-frequZLyr[count.fields(textConnection(frequZLyr$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZLyr,"CnvregionswithmeanarenosaZapaLyrata_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZLyr<-read.xlsx("CnvregionswithmeanarenosaZapaLyrata_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigZLyr[,2]<-as.numeric(as.character(frequbigZLyr[,2]))
frequbigZLyr[,3]<-as.numeric(as.character(frequbigZLyr[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZLyr_GRange<-GRanges(seqnames=tolower(frequbigZLyr$Scaffold),ranges=IRanges(start=frequbigZLyr$start,end=frequbigZLyr$end))
values(frequbigZLyr_GRange)<-frequbigZLyr[,6:11]

testZLyr<-findOverlaps(frequbigZLyr_GRange,lyrGenes)
test2ZLyr<-pintersect(frequbigZLyr_GRange[queryHits(testZLyr)],lyrGenes[subjectHits(testZLyr)])
test3ZLyr<-as.data.frame(test2ZLyr)
frequbigZLyr_lyrGenes=mergeByOverlaps(frequbigZLyr_GRange,lyrGenes,type=c("any"))
frequbigZLyr_lyrGenes_df=data.frame(as.data.frame(frequbigZLyr_lyrGenes$frequbigZLyr_GRange),as.data.frame(frequbigZLyr_lyrGenes$lyrGenes),test3ZLyr$width)
colnames(frequbigZLyr_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZLyr_lyrGenes_df<-frequbigZLyr_lyrGenes_df[(frequbigZLyr_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZLyr_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZLyr_lyrGenes_df,"ZLyr_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZLyr_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZLyr_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZLyr_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup=merge(frequbigZLyr_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZLyr.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZLyr.xlsx",overwrite=T)

Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZLyr_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaZapaLyratacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Genesarenosacnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

LZar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LZar.xlsx",1)
ZLyr<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZLyr.xlsx",1)

LZar_unique<-LZar[!duplicated(LZar$Lyr_Gene,LZar$Copy_start,LZar$Copy_end),]
ZLyr_unique<-ZLyr[!duplicated(ZLyr$Lyr_Gene,ZLyr$Copy_start,ZLyr$Copy_end),]

LZar_unique$Lyr_Gene<-as.character(LZar_unique$Lyr_Gene)
ZLyr_unique$Lyr_Gene<-as.character(ZLyr_unique$Lyr_Gene)
LZar_unique$CN_class<-as.character(LZar_unique$CN_class)
ZLyr_unique$CN_class<-as.character(ZLyr_unique$CN_class)

LZaronly<-LZar_unique[!((LZar_unique$Lyr_Gene%in%ZLyr_unique$Lyr_Gene)&(LZar_unique$CN_class==ZLyr_unique$CN_class[match(LZar_unique$Lyr_Gene,ZLyr_unique$Lyr_Gene)])),]
ZLyronly<-ZLyr_unique[!((ZLyr_unique$Lyr_Gene%in%LZar_unique$Lyr_Gene)&(ZLyr_unique$CN_class==LZar_unique$CN_class[match(ZLyr_unique$Lyr_Gene,LZar_unique$Lyr_Gene)])),]

write.table(LZaronly,"LZar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(ZLyronly,"ZLyr_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











