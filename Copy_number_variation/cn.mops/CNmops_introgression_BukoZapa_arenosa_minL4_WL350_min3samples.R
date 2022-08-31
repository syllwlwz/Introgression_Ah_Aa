#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesBukoarenosa <- scan("Buko_arenosa_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesBukoarenosa)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Bukoarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesBukoarenosa,sampleNames=c("Buko_a21","Buko_a22","Buko_a28","Buko_a29","Buko_a30","Buko_a31","Buko_a34b","Buko_a34"),refSeqName=sort(scaffolds),WL=350) 
Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapaarenosa,sampleNames=c("Zapa002a01","Zapa002a02","Zapa002a03","Zapa002a04","Zapa002a05","Zapa002a07","Zapa004a11","Zapa008a19","Zapa_11a03x12a01_l","Zapa_a09xa03_b"),refSeqName=sort(scaffolds),WL=350)
BZar <- referencecn.mops(cases=Bukoarenosa,controls=Zapaarenosa, minWidth = 4, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSBZar <- calcIntegerCopyNumbers(BZar)

ZBar <- referencecn.mops(cases=Zapaarenosa,controls=Bukoarenosa, minWidth = 4, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSZBar <- calcIntegerCopyNumbers(ZBar)


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
regionsBZar<-cnvr(resCNMOPSBZar)
singleBZar<-cnvs(resCNMOPSBZar)
testBZar<-findOverlaps(regionsBZar,singleBZar)
test2BZar<- DataFrame(splitAsList(singleBZar$sampleName[subjectHits(testBZar)], queryHits(testBZar)),splitAsList(singleBZar$median[subjectHits(testBZar)], queryHits(testBZar)),splitAsList(singleBZar$mean[subjectHits(testBZar)], queryHits(testBZar)),splitAsList(singleBZar$CN[subjectHits(testBZar)], queryHits(testBZar)))
colnames(test2BZar)<-c("sampleNames","mean","median","CN")
test7BZar<-apply(as.data.frame(test2BZar),1,dupl_cov)
mcols(regionsBZar)<-as.data.frame(as.matrix(test7BZar))

dataBZar<-as.data.frame(regionsBZar)

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

adddataBZar<-do.call("rbind",apply(dataBZar,1,data_change))
dataBZar<-dataBZar[,-6]
frequBZar<-cbind(dataBZar,adddataBZar)

names(frequBZar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3BZar<-frequBZar[count.fields(textConnection(frequBZar$sampleNames),sep=",")>=5,]
frequ4BZar<-frequBZar[count.fields(textConnection(frequBZar$sampleNames),sep=",")>=3,]

require(openxlsx)
write.xlsx(frequ3BZar,"CnvregionswithmeanarenosaBukoZapaarenosa_minRead6_WL350_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
write.xlsx(frequ4BZar,"CnvregionswithmeanarenosaBukoZapaarenosa_minRead6_WL350_minL4_min3samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigBZar<-read.xlsx("CnvregionswithmeanarenosaBukoZapaarenosa_minRead6_WL350_minL4_min3samples.xlsx",1)
frequbigBZar[,2]<-as.numeric(as.character(frequbigBZar[,2]))
frequbigBZar[,3]<-as.numeric(as.character(frequbigBZar[,3]))

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
frequbigBZar_GRange<-GRanges(seqnames=tolower(frequbigBZar$Scaffold),ranges=IRanges(start=frequbigBZar$start,end=frequbigBZar$end))
values(frequbigBZar_GRange)<-frequbigBZar[,6:11]

testBZar<-findOverlaps(frequbigBZar_GRange,lyrGenes)
test2BZar<-pintersect(frequbigBZar_GRange[queryHits(testBZar)],lyrGenes[subjectHits(testBZar)])
test3BZar<-as.data.frame(test2BZar)
frequbigBZar_lyrGenes=mergeByOverlaps(frequbigBZar_GRange,lyrGenes,type=c("any"))
frequbigBZar_lyrGenes_df=data.frame(as.data.frame(frequbigBZar_lyrGenes$frequbigBZar_GRange),as.data.frame(frequbigBZar_lyrGenes$lyrGenes),test3BZar$width)
colnames(frequbigBZar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigBZar_lyrGenes_df<-frequbigBZar_lyrGenes_df[(frequbigBZar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigBZar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigBZar_lyrGenes_df,"BZar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigBZar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigBZar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigBZar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup=merge(frequbigBZar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BZar_min3samples.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BZar_min3samples.xlsx",overwrite=T)

Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigBZar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaBukoZapaarenosacnvregions_minRead6_WL350_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast3samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Buko arenosa#
#################################################

regionsZBar<-cnvr(resCNMOPSZBar)
singleZBar<-cnvs(resCNMOPSZBar)
testZBar<-findOverlaps(regionsZBar,singleZBar)
test2ZBar<- DataFrame(splitAsList(singleZBar$sampleName[subjectHits(testZBar)], queryHits(testZBar)),splitAsList(singleZBar$median[subjectHits(testZBar)], queryHits(testZBar)),splitAsList(singleZBar$mean[subjectHits(testZBar)], queryHits(testZBar)),splitAsList(singleZBar$CN[subjectHits(testZBar)], queryHits(testZBar)))
colnames(test2ZBar)<-c("sampleNames","mean","median","CN")
test7ZBar<-apply(as.data.frame(test2ZBar),1,dupl_cov)
mcols(regionsZBar)<-as.data.frame(as.matrix(test7ZBar))

dataZBar<-as.data.frame(regionsZBar)

adddataZBar<-do.call("rbind",apply(dataZBar,1,data_change))
dataZBar<-dataZBar[,-6]
frequZBar<-cbind(dataZBar,adddataZBar)

names(frequZBar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZBar<-frequZBar[count.fields(textConnection(frequZBar$sampleNames),sep=",")>=5,]
frequ4ZBar<-frequZBar[count.fields(textConnection(frequZBar$sampleNames),sep=",")>=4,]
require(openxlsx)
write.xlsx(frequ3ZBar,"CnvregionswithmeanarenosaZapaBukoarenosa_minRead6_WL350_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
write.xlsx(frequ4ZBar,"CnvregionswithmeanarenosaZapaBukoarenosa_minRead6_WL350_minL4_min3samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZBar<-read.xlsx("CnvregionswithmeanarenosaZapaBukoarenosa_minRead6_WL350_minL4_min3samples.xlsx",1)
frequbigZBar[,2]<-as.numeric(as.character(frequbigZBar[,2]))
frequbigZBar[,3]<-as.numeric(as.character(frequbigZBar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZBar_GRange<-GRanges(seqnames=tolower(frequbigZBar$Scaffold),ranges=IRanges(start=frequbigZBar$start,end=frequbigZBar$end))
values(frequbigZBar_GRange)<-frequbigZBar[,6:11]

testZBar<-findOverlaps(frequbigZBar_GRange,lyrGenes)
test2ZBar<-pintersect(frequbigZBar_GRange[queryHits(testZBar)],lyrGenes[subjectHits(testZBar)])
test3ZBar<-as.data.frame(test2ZBar)
frequbigZBar_lyrGenes=mergeByOverlaps(frequbigZBar_GRange,lyrGenes,type=c("any"))
frequbigZBar_lyrGenes_df=data.frame(as.data.frame(frequbigZBar_lyrGenes$frequbigZBar_GRange),as.data.frame(frequbigZBar_lyrGenes$lyrGenes),test3ZBar$width)
colnames(frequbigZBar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZBar_lyrGenes_df<-frequbigZBar_lyrGenes_df[(frequbigZBar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZBar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZBar_lyrGenes_df,"ZBar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZBar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZBar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZBar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup=merge(frequbigZBar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZBar_min3samples.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZBar_min3samples.xlsx",overwrite=T)

Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZBar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaZapaBukoarenosacnvregions_minRead6_WL350_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast3samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

BZar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_BZar_min3samples.xlsx",1)
ZBar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZBar_min3samples.xlsx",1)

BZar_unique<-BZar[!duplicated(BZar$Lyr_Gene,BZar$Copy_start,BZar$Copy_end),]
ZBar_unique<-ZBar[!duplicated(ZBar$Lyr_Gene,ZBar$Copy_start,ZBar$Copy_end),]

BZar_unique$Lyr_Gene<-as.character(BZar_unique$Lyr_Gene)
ZBar_unique$Lyr_Gene<-as.character(ZBar_unique$Lyr_Gene)
BZar_unique$CN_class<-as.character(BZar_unique$CN_class)
ZBar_unique$CN_class<-as.character(ZBar_unique$CN_class)

BZaronly<-BZar_unique[!((BZar_unique$Lyr_Gene%in%ZBar_unique$Lyr_Gene)&(BZar_unique$CN_class==ZBar_unique$CN_class[match(BZar_unique$Lyr_Gene,ZBar_unique$Lyr_Gene)])),]
ZBaronly<-ZBar_unique[!((ZBar_unique$Lyr_Gene%in%BZar_unique$Lyr_Gene)&(ZBar_unique$CN_class==BZar_unique$CN_class[match(ZBar_unique$Lyr_Gene,BZar_unique$Lyr_Gene)])),]

write.table(BZaronly,"BZar_cnvs_only_strict_minL4_WL350_min3samples.table",row.names=F,sep="\t")
write.table(ZBaronly,"ZBar_cnvs_only_strict_minL4_WL350_min3samples.table",row.names=F,sep="\t")











