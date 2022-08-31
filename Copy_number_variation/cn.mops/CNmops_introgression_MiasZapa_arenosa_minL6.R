#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesMiasarenosa <- scan("Mias_arenosa_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesMiasarenosa)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Miasarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesMiasarenosa,sampleNames=c("Mias001a18","Mias003a09","Mias003a10","Mias003a11","Mias003a13","Mias003a15","Mias003a16","Mias_19a03x21a01_d","Mias_23a03x22a03_k","Mias_a12xa05_l","Mias_a42","Mias_a43","Mias_a44","Mias_a45","Mias_a46","Mias_a47","Mias_a48","Mias_a50","Mias_a55"),refSeqName=sort(scaffolds),WL=500) 
Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapaarenosa,sampleNames=c("Zapa002a01","Zapa002a02","Zapa002a03","Zapa002a04","Zapa002a05","Zapa002a07","Zapa004a11","Zapa008a19","Zapa_11a03x12a01_l","Zapa_a09xa03_b"),refSeqName=sort(scaffolds),WL=500) 
MZar <- referencecn.mops(cases=Miasarenosa,controls=Zapaarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSMZar <- calcIntegerCopyNumbers(MZar)

ZMar <- referencecn.mops(cases=Zapaarenosa,controls=Miasarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSZMar <- calcIntegerCopyNumbers(ZMar)


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
regionsMZar<-cnvr(resCNMOPSMZar)
singleMZar<-cnvs(resCNMOPSMZar)
testMZar<-findOverlaps(regionsMZar,singleMZar)
test2MZar<- DataFrame(splitAsList(singleMZar$sampleName[subjectHits(testMZar)], queryHits(testMZar)),splitAsList(singleMZar$median[subjectHits(testMZar)], queryHits(testMZar)),splitAsList(singleMZar$mean[subjectHits(testMZar)], queryHits(testMZar)),splitAsList(singleMZar$CN[subjectHits(testMZar)], queryHits(testMZar)))
colnames(test2MZar)<-c("sampleNames","mean","median","CN")
test7MZar<-apply(as.data.frame(test2MZar),1,dupl_cov)
mcols(regionsMZar)<-as.data.frame(as.matrix(test7MZar))

dataMZar<-as.data.frame(regionsMZar)

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

adddataMZar<-do.call("rbind",apply(dataMZar,1,data_change))
dataMZar<-dataMZar[,-6]
frequMZar<-cbind(dataMZar,adddataMZar)

names(frequMZar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3MZar<-frequMZar[count.fields(textConnection(frequMZar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3MZar,"CnvregionswithmeanarenosaMiasZapaarenosa_minRead6_WL500_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigMZar<-read.xlsx("CnvregionswithmeanarenosaMiasZapaarenosa_minRead6_WL500_minL6_min5samples.xlsx",1)
frequbigMZar[,2]<-as.numeric(as.character(frequbigMZar[,2]))
frequbigMZar[,3]<-as.numeric(as.character(frequbigMZar[,3]))

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
frequbigMZar_GRange<-GRanges(seqnames=tolower(frequbigMZar$Scaffold),ranges=IRanges(start=frequbigMZar$start,end=frequbigMZar$end))
values(frequbigMZar_GRange)<-frequbigMZar[,6:11]

testMZar<-findOverlaps(frequbigMZar_GRange,lyrGenes)
test2MZar<-pintersect(frequbigMZar_GRange[queryHits(testMZar)],lyrGenes[subjectHits(testMZar)])
test3MZar<-as.data.frame(test2MZar)
frequbigMZar_lyrGenes=mergeByOverlaps(frequbigMZar_GRange,lyrGenes,type=c("any"))
frequbigMZar_lyrGenes_df=data.frame(as.data.frame(frequbigMZar_lyrGenes$frequbigMZar_GRange),as.data.frame(frequbigMZar_lyrGenes$lyrGenes),test3MZar$width)
colnames(frequbigMZar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigMZar_lyrGenes_df<-frequbigMZar_lyrGenes_df[(frequbigMZar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigMZar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigMZar_lyrGenes_df,"MZar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigMZar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigMZar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigMZar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup=merge(frequbigMZar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MZar_minL6.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MZar_minL6.xlsx",overwrite=T)

Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigMZar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaMiasZapaarenosacnvregions_minRead6_WL500_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Mias arenosa#
#################################################

regionsZMar<-cnvr(resCNMOPSZMar)
singleZMar<-cnvs(resCNMOPSZMar)
testZMar<-findOverlaps(regionsZMar,singleZMar)
test2ZMar<- DataFrame(splitAsList(singleZMar$sampleName[subjectHits(testZMar)], queryHits(testZMar)),splitAsList(singleZMar$median[subjectHits(testZMar)], queryHits(testZMar)),splitAsList(singleZMar$mean[subjectHits(testZMar)], queryHits(testZMar)),splitAsList(singleZMar$CN[subjectHits(testZMar)], queryHits(testZMar)))
colnames(test2ZMar)<-c("sampleNames","mean","median","CN")
test7ZMar<-apply(as.data.frame(test2ZMar),1,dupl_cov)
mcols(regionsZMar)<-as.data.frame(as.matrix(test7ZMar))

dataZMar<-as.data.frame(regionsZMar)

adddataZMar<-do.call("rbind",apply(dataZMar,1,data_change))
dataZMar<-dataZMar[,-6]
frequZMar<-cbind(dataZMar,adddataZMar)

names(frequZMar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZMar<-frequZMar[count.fields(textConnection(frequZMar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZMar,"CnvregionswithmeanarenosaZapaMiasarenosa_minRead6_WL500_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZMar<-read.xlsx("CnvregionswithmeanarenosaZapaMiasarenosa_minRead6_WL500_minL6_min5samples.xlsx",1)
frequbigZMar[,2]<-as.numeric(as.character(frequbigZMar[,2]))
frequbigZMar[,3]<-as.numeric(as.character(frequbigZMar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZMar_GRange<-GRanges(seqnames=tolower(frequbigZMar$Scaffold),ranges=IRanges(start=frequbigZMar$start,end=frequbigZMar$end))
values(frequbigZMar_GRange)<-frequbigZMar[,6:11]

testZMar<-findOverlaps(frequbigZMar_GRange,lyrGenes)
test2ZMar<-pintersect(frequbigZMar_GRange[queryHits(testZMar)],lyrGenes[subjectHits(testZMar)])
test3ZMar<-as.data.frame(test2ZMar)
frequbigZMar_lyrGenes=mergeByOverlaps(frequbigZMar_GRange,lyrGenes,type=c("any"))
frequbigZMar_lyrGenes_df=data.frame(as.data.frame(frequbigZMar_lyrGenes$frequbigZMar_GRange),as.data.frame(frequbigZMar_lyrGenes$lyrGenes),test3ZMar$width)
colnames(frequbigZMar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZMar_lyrGenes_df<-frequbigZMar_lyrGenes_df[(frequbigZMar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZMar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZMar_lyrGenes_df,"ZMar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZMar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZMar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZMar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup=merge(frequbigZMar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZMar_minL6.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZMar_minL6.xlsx",overwrite=T)

Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZMar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaZapaMiasarenosacnvregions_minRead6_WL500_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

MZar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MZar_minL6.xlsx",1)
ZMar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZMar_minL6.xlsx",1)

MZar_unique<-MZar[!duplicated(MZar$Lyr_Gene,MZar$Copy_start,MZar$Copy_end),]
ZMar_unique<-ZMar[!duplicated(ZMar$Lyr_Gene,ZMar$Copy_start,ZMar$Copy_end),]

MZar_unique$Lyr_Gene<-as.character(MZar_unique$Lyr_Gene)
ZMar_unique$Lyr_Gene<-as.character(ZMar_unique$Lyr_Gene)
MZar_unique$CN_class<-as.character(MZar_unique$CN_class)
ZMar_unique$CN_class<-as.character(ZMar_unique$CN_class)

MZaronly<-MZar_unique[!((MZar_unique$Lyr_Gene%in%ZMar_unique$Lyr_Gene)&(MZar_unique$CN_class==ZMar_unique$CN_class[match(MZar_unique$Lyr_Gene,ZMar_unique$Lyr_Gene)])),]
ZMaronly<-ZMar_unique[!((ZMar_unique$Lyr_Gene%in%MZar_unique$Lyr_Gene)&(ZMar_unique$CN_class==MZar_unique$CN_class[match(ZMar_unique$Lyr_Gene,MZar_unique$Lyr_Gene)])),]

write.table(MZaronly,"MZar_cnvs_only_strict_minL6.table",row.names=F,sep="\t")
write.table(ZMaronly,"ZMar_cnvs_only_strict_minL6.table",row.names=F,sep="\t")











