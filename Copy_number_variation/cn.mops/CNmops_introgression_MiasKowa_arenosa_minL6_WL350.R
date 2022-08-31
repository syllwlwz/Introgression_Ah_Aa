#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesMiasarenosa <- scan("Mias_arenosa_full.list",character(),quote="")
BAMFilesKowaarenosa <- scan("Kowa_arenosa_full.list",character(),quote="")

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

Miasarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesMiasarenosa,sampleNames=c("Mias001a18","Mias003a09","Mias003a10","Mias003a11","Mias003a13","Mias003a15","Mias003a16","Mias_19a03x21a01_d","Mias_23a03x22a03_k","Mias_a12xa05_l","Mias_a42","Mias_a43","Mias_a44","Mias_a45","Mias_a46","Mias_a47","Mias_a48","Mias_a50","Mias_a55"),refSeqName=sort(scaffolds),WL=350) 
Kowaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesKowaarenosa,sampleNames=c("Kowa001a04","Kowa001a05","Kowa001a06","Kowa001a07","Kowa001a08","Kowa001a09","Kowa001a11","Kowa001a12"),refSeqName=sort(scaffolds),WL=350) 
MKar <- referencecn.mops(cases=Miasarenosa,controls=Kowaarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSMKar <- calcIntegerCopyNumbers(MKar)

KMar <- referencecn.mops(cases=Kowaarenosa,controls=Miasarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSKMar <- calcIntegerCopyNumbers(KMar)


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
regionsMKar<-cnvr(resCNMOPSMKar)
singleMKar<-cnvs(resCNMOPSMKar)
testMKar<-findOverlaps(regionsMKar,singleMKar)
test2MKar<- DataFrame(splitAsList(singleMKar$sampleName[subjectHits(testMKar)], queryHits(testMKar)),splitAsList(singleMKar$median[subjectHits(testMKar)], queryHits(testMKar)),splitAsList(singleMKar$mean[subjectHits(testMKar)], queryHits(testMKar)),splitAsList(singleMKar$CN[subjectHits(testMKar)], queryHits(testMKar)))
colnames(test2MKar)<-c("sampleNames","mean","median","CN")
test7MKar<-apply(as.data.frame(test2MKar),1,dupl_cov)
mcols(regionsMKar)<-as.data.frame(as.matrix(test7MKar))

dataMKar<-as.data.frame(regionsMKar)

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

adddataMKar<-do.call("rbind",apply(dataMKar,1,data_change))
dataMKar<-dataMKar[,-6]
frequMKar<-cbind(dataMKar,adddataMKar)

names(frequMKar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3MKar<-frequMKar[count.fields(textConnection(frequMKar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3MKar,"CnvregionswithmeanarenosaMiasKowaarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigMKar<-read.xlsx("CnvregionswithmeanarenosaMiasKowaarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigMKar[,2]<-as.numeric(as.character(frequbigMKar[,2]))
frequbigMKar[,3]<-as.numeric(as.character(frequbigMKar[,3]))

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
frequbigMKar_GRange<-GRanges(seqnames=tolower(frequbigMKar$Scaffold),ranges=IRanges(start=frequbigMKar$start,end=frequbigMKar$end))
values(frequbigMKar_GRange)<-frequbigMKar[,6:11]

testMKar<-findOverlaps(frequbigMKar_GRange,lyrGenes)
test2MKar<-pintersect(frequbigMKar_GRange[queryHits(testMKar)],lyrGenes[subjectHits(testMKar)])
test3MKar<-as.data.frame(test2MKar)
frequbigMKar_lyrGenes=mergeByOverlaps(frequbigMKar_GRange,lyrGenes,type=c("any"))
frequbigMKar_lyrGenes_df=data.frame(as.data.frame(frequbigMKar_lyrGenes$frequbigMKar_GRange),as.data.frame(frequbigMKar_lyrGenes$lyrGenes),test3MKar$width)
colnames(frequbigMKar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigMKar_lyrGenes_df<-frequbigMKar_lyrGenes_df[(frequbigMKar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigMKar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigMKar_lyrGenes_df,"MKar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx",1)

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigMKar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigMKar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Alyr_ID",all.x=TRUE)
colnames(frequbigMKar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup=merge(frequbigMKar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MKar_minL6_WL350.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MKar_minL6_WL350.xlsx",overwrite=T)

Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigMKar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaMiasKowaarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Mias arenosa#
#################################################

regionsKMar<-cnvr(resCNMOPSKMar)
singleKMar<-cnvs(resCNMOPSKMar)
testKMar<-findOverlaps(regionsKMar,singleKMar)
test2KMar<- DataFrame(splitAsList(singleKMar$sampleName[subjectHits(testKMar)], queryHits(testKMar)),splitAsList(singleKMar$median[subjectHits(testKMar)], queryHits(testKMar)),splitAsList(singleKMar$mean[subjectHits(testKMar)], queryHits(testKMar)),splitAsList(singleKMar$CN[subjectHits(testKMar)], queryHits(testKMar)))
colnames(test2KMar)<-c("sampleNames","mean","median","CN")
test7KMar<-apply(as.data.frame(test2KMar),1,dupl_cov)
mcols(regionsKMar)<-as.data.frame(as.matrix(test7KMar))

dataKMar<-as.data.frame(regionsKMar)

adddataKMar<-do.call("rbind",apply(dataKMar,1,data_change))
dataKMar<-dataKMar[,-6]
frequKMar<-cbind(dataKMar,adddataKMar)

names(frequKMar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KMar<-frequKMar[count.fields(textConnection(frequKMar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KMar,"CnvregionswithmeanarenosaKowaMiasarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKMar<-read.xlsx("CnvregionswithmeanarenosaKowaMiasarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigKMar[,2]<-as.numeric(as.character(frequbigKMar[,2]))
frequbigKMar[,3]<-as.numeric(as.character(frequbigKMar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKMar_GRange<-GRanges(seqnames=tolower(frequbigKMar$Scaffold),ranges=IRanges(start=frequbigKMar$start,end=frequbigKMar$end))
values(frequbigKMar_GRange)<-frequbigKMar[,6:11]

testKMar<-findOverlaps(frequbigKMar_GRange,lyrGenes)
test2KMar<-pintersect(frequbigKMar_GRange[queryHits(testKMar)],lyrGenes[subjectHits(testKMar)])
test3KMar<-as.data.frame(test2KMar)
frequbigKMar_lyrGenes=mergeByOverlaps(frequbigKMar_GRange,lyrGenes,type=c("any"))
frequbigKMar_lyrGenes_df=data.frame(as.data.frame(frequbigKMar_lyrGenes$frequbigKMar_GRange),as.data.frame(frequbigKMar_lyrGenes$lyrGenes),test3KMar$width)
colnames(frequbigKMar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKMar_lyrGenes_df<-frequbigKMar_lyrGenes_df[(frequbigKMar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKMar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKMar_lyrGenes_df,"KMar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigKMar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKMar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Alyr_ID",all.x=TRUE)
colnames(frequbigKMar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup=merge(frequbigKMar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KMar_minL6_WL350.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KMar_minL6_WL350.xlsx",overwrite=T)

Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKMar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaKowaMiasarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

MKar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MKar_minL6_WL350.xlsx",1)
KMar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KMar_minL6_WL350.xlsx",1)

MKar_unique<-MKar[!duplicated(MKar$Lyr_Gene,MKar$Copy_start,MKar$Copy_end),]
KMar_unique<-KMar[!duplicated(KMar$Lyr_Gene,KMar$Copy_start,KMar$Copy_end),]

MKar_unique$Lyr_Gene<-as.character(MKar_unique$Lyr_Gene)
KMar_unique$Lyr_Gene<-as.character(KMar_unique$Lyr_Gene)
MKar_unique$CN_class<-as.character(MKar_unique$CN_class)
KMar_unique$CN_class<-as.character(KMar_unique$CN_class)

MKaronly<-MKar_unique[!((MKar_unique$Lyr_Gene%in%KMar_unique$Lyr_Gene)&(MKar_unique$CN_class==KMar_unique$CN_class[match(MKar_unique$Lyr_Gene,KMar_unique$Lyr_Gene)])),]
KMaronly<-KMar_unique[!((KMar_unique$Lyr_Gene%in%MKar_unique$Lyr_Gene)&(KMar_unique$CN_class==MKar_unique$CN_class[match(KMar_unique$Lyr_Gene,MKar_unique$Lyr_Gene)])),]

write.table(MKaronly,"MKar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(KMaronly,"KMar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











