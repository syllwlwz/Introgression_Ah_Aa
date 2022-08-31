#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesBukoarenosa <- scan("Buko_arenosa_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesBukoarenosa)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Bukoarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesBukoarenosa,sampleNames=c("Buko_a21","Buko_a22","Buko_a28","Buko_a29","Buko_a30","Buko_a31","Buko_a34b","Buko_a34"),refSeqName=sort(scaffolds),WL=350) 
Kowaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesKowaarenosa,sampleNames=c("Kowa001a04","Kowa001a05","Kowa001a06","Kowa001a07","Kowa001a08","Kowa001a09","Kowa001a11","Kowa001a12"),refSeqName=sort(scaffolds),WL=350) 
BKar <- referencecn.mops(cases=Bukoarenosa,controls=Kowaarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSBKar <- calcIntegerCopyNumbers(BKar)

KBar <- referencecn.mops(cases=Kowaarenosa,controls=Bukoarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSKBar <- calcIntegerCopyNumbers(KBar)


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
regionsBKar<-cnvr(resCNMOPSBKar)
singleBKar<-cnvs(resCNMOPSBKar)
testBKar<-findOverlaps(regionsBKar,singleBKar)
test2BKar<- DataFrame(splitAsList(singleBKar$sampleName[subjectHits(testBKar)], queryHits(testBKar)),splitAsList(singleBKar$median[subjectHits(testBKar)], queryHits(testBKar)),splitAsList(singleBKar$mean[subjectHits(testBKar)], queryHits(testBKar)),splitAsList(singleBKar$CN[subjectHits(testBKar)], queryHits(testBKar)))
colnames(test2BKar)<-c("sampleNames","mean","median","CN")
test7BKar<-apply(as.data.frame(test2BKar),1,dupl_cov)
mcols(regionsBKar)<-as.data.frame(as.matrix(test7BKar))

dataBKar<-as.data.frame(regionsBKar)

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

adddataBKar<-do.call("rbind",apply(dataBKar,1,data_change))
dataBKar<-dataBKar[,-6]
frequBKar<-cbind(dataBKar,adddataBKar)

names(frequBKar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3BKar<-frequBKar[count.fields(textConnection(frequBKar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3BKar,"CnvregionswithmeanarenosaBukoKowaarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
write.xlsx(frequBKar,"CnvregionswithmeanarenosaBukoKowaarenosa_minRead6_WL350_minL6_test.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigBKar<-read.xlsx("CnvregionswithmeanarenosaBukoKowaarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigBKar[,2]<-as.numeric(as.character(frequbigBKar[,2]))
frequbigBKar[,3]<-as.numeric(as.character(frequbigBKar[,3]))

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
frequbigBKar_GRange<-GRanges(seqnames=tolower(frequbigBKar$Scaffold),ranges=IRanges(start=frequbigBKar$start,end=frequbigBKar$end))
values(frequbigBKar_GRange)<-frequbigBKar[,6:11]

testBKar<-findOverlaps(frequbigBKar_GRange,lyrGenes)
test2BKar<-pintersect(frequbigBKar_GRange[queryHits(testBKar)],lyrGenes[subjectHits(testBKar)])
test3BKar<-as.data.frame(test2BKar)
frequbigBKar_lyrGenes=mergeByOverlaps(frequbigBKar_GRange,lyrGenes,type=c("any"))
frequbigBKar_lyrGenes_df=data.frame(as.data.frame(frequbigBKar_lyrGenes$frequbigBKar_GRange),as.data.frame(frequbigBKar_lyrGenes$lyrGenes),test3BKar$width)
colnames(frequbigBKar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigBKar_lyrGenes_df<-frequbigBKar_lyrGenes_df[(frequbigBKar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigBKar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigBKar_lyrGenes_df,"BKar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigBKar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigBKar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigBKar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup=merge(frequbigBKar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BKar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BKar.xlsx",overwrite=T)

Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigBKar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaBukoKowaarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Genesarenosacnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Buko arenosa#
#################################################

regionsKBar<-cnvr(resCNMOPSKBar)
singleKBar<-cnvs(resCNMOPSKBar)
testKBar<-findOverlaps(regionsKBar,singleKBar)
test2KBar<- DataFrame(splitAsList(singleKBar$sampleName[subjectHits(testKBar)], queryHits(testKBar)),splitAsList(singleKBar$median[subjectHits(testKBar)], queryHits(testKBar)),splitAsList(singleKBar$mean[subjectHits(testKBar)], queryHits(testKBar)),splitAsList(singleKBar$CN[subjectHits(testKBar)], queryHits(testKBar)))
colnames(test2KBar)<-c("sampleNames","mean","median","CN")
test7KBar<-apply(as.data.frame(test2KBar),1,dupl_cov)
mcols(regionsKBar)<-as.data.frame(as.matrix(test7KBar))

dataKBar<-as.data.frame(regionsKBar)

adddataKBar<-do.call("rbind",apply(dataKBar,1,data_change))
dataKBar<-dataKBar[,-6]
frequKBar<-cbind(dataKBar,adddataKBar)

names(frequKBar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KBar<-frequKBar[count.fields(textConnection(frequKBar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KBar,"CnvregionswithmeanarenosaKowaBukoarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKBar<-read.xlsx("CnvregionswithmeanarenosaKowaBukoarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigKBar[,2]<-as.numeric(as.character(frequbigKBar[,2]))
frequbigKBar[,3]<-as.numeric(as.character(frequbigKBar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKBar_GRange<-GRanges(seqnames=tolower(frequbigKBar$Scaffold),ranges=IRanges(start=frequbigKBar$start,end=frequbigKBar$end))
values(frequbigKBar_GRange)<-frequbigKBar[,6:11]

testKBar<-findOverlaps(frequbigKBar_GRange,lyrGenes)
test2KBar<-pintersect(frequbigKBar_GRange[queryHits(testKBar)],lyrGenes[subjectHits(testKBar)])
test3KBar<-as.data.frame(test2KBar)
frequbigKBar_lyrGenes=mergeByOverlaps(frequbigKBar_GRange,lyrGenes,type=c("any"))
frequbigKBar_lyrGenes_df=data.frame(as.data.frame(frequbigKBar_lyrGenes$frequbigKBar_GRange),as.data.frame(frequbigKBar_lyrGenes$lyrGenes),test3KBar$width)
colnames(frequbigKBar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKBar_lyrGenes_df<-frequbigKBar_lyrGenes_df[(frequbigKBar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKBar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKBar_lyrGenes_df,"KBar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigKBar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKBar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKBar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup=merge(frequbigKBar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KBar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KBar.xlsx",overwrite=T)

Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKBar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaKowaBukoarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Genesarenosacnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

BKar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BKar.xlsx",1)
KBar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KBar.xlsx",1)

BKar_unique<-BKar[!duplicated(BKar$Lyr_Gene,BKar$Copy_start,BKar$Copy_end),]
KBar_unique<-KBar[!duplicated(KBar$Lyr_Gene,KBar$Copy_start,KBar$Copy_end),]

BKar_unique$Lyr_Gene<-as.character(BKar_unique$Lyr_Gene)
KBar_unique$Lyr_Gene<-as.character(KBar_unique$Lyr_Gene)
BKar_unique$CN_class<-as.character(BKar_unique$CN_class)
KBar_unique$CN_class<-as.character(KBar_unique$CN_class)

BKaronly<-BKar_unique[!((BKar_unique$Lyr_Gene%in%KBar_unique$Lyr_Gene)&(BKar_unique$CN_class==KBar_unique$CN_class[match(BKar_unique$Lyr_Gene,KBar_unique$Lyr_Gene)])),]
KBaronly<-KBar_unique[!((KBar_unique$Lyr_Gene%in%BKar_unique$Lyr_Gene)&(KBar_unique$CN_class==BKar_unique$CN_class[match(KBar_unique$Lyr_Gene,BKar_unique$Lyr_Gene)])),]

write.table(BKaronly,"BKar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(KBaronly,"KBar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











