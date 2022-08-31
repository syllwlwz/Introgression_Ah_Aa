#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesPiekarenosa <- scan("Piek_arenosa_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesPiekarenosa)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Piekarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesPiekarenosa,sampleNames=c("Piek_a10","Piek_a12","Piek_a13","Piek_a14","Piek_a15","Piek_a16","Piek_a17","Piek_a18","Piek_a19","Piek_a1","Piek_a20-1","Piek_a20-2","Piek_a22","Piek_a2","Piek_a4","Piek_a5","Piek_a7","Piek_a8","Piek_a9","Piek_h14_1","Piek_h6_2"),refSeqName=sort(scaffolds),WL=300) 
Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapaarenosa,sampleNames=c("Zapa002a01","Zapa002a02","Zapa002a03","Zapa002a04","Zapa002a05","Zapa002a07","Zapa004a11","Zapa008a19","Zapa_11a03x12a01_l","Zapa_a09xa03_b"),refSeqName=sort(scaffolds),WL=300) 
PZar <- referencecn.mops(cases=Piekarenosa,controls=Zapaarenosa, minWidth = 4, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSPZar <- calcIntegerCopyNumbers(PZar)

ZPar <- referencecn.mops(cases=Zapaarenosa,controls=Piekarenosa, minWidth = 4, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSZPar <- calcIntegerCopyNumbers(ZPar)


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
regionsPZar<-cnvr(resCNMOPSPZar)
singlePZar<-cnvs(resCNMOPSPZar)
testPZar<-findOverlaps(regionsPZar,singlePZar)
test2PZar<- DataFrame(splitAsList(singlePZar$sampleName[subjectHits(testPZar)], queryHits(testPZar)),splitAsList(singlePZar$median[subjectHits(testPZar)], queryHits(testPZar)),splitAsList(singlePZar$mean[subjectHits(testPZar)], queryHits(testPZar)),splitAsList(singlePZar$CN[subjectHits(testPZar)], queryHits(testPZar)))
colnames(test2PZar)<-c("sampleNames","mean","median","CN")
test7PZar<-apply(as.data.frame(test2PZar),1,dupl_cov)
mcols(regionsPZar)<-as.data.frame(as.matrix(test7PZar))

dataPZar<-as.data.frame(regionsPZar)

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

adddataPZar<-do.call("rbind",apply(dataPZar,1,data_change))
dataPZar<-dataPZar[,-6]
frequPZar<-cbind(dataPZar,adddataPZar)

names(frequPZar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3PZar<-frequPZar[count.fields(textConnection(frequPZar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3PZar,"CnvregionswithmeanarenosaPiekZapaarenosa_minRead6_WL300_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigPZar<-read.xlsx("CnvregionswithmeanarenosaPiekZapaarenosa_minRead6_WL300_minL4_min5samples.xlsx",1)
frequbigPZar[,2]<-as.numeric(as.character(frequbigPZar[,2]))
frequbigPZar[,3]<-as.numeric(as.character(frequbigPZar[,3]))

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
frequbigPZar_GRange<-GRanges(seqnames=tolower(frequbigPZar$Scaffold),ranges=IRanges(start=frequbigPZar$start,end=frequbigPZar$end))
values(frequbigPZar_GRange)<-frequbigPZar[,6:11]

testPZar<-findOverlaps(frequbigPZar_GRange,lyrGenes)
test2PZar<-pintersect(frequbigPZar_GRange[queryHits(testPZar)],lyrGenes[subjectHits(testPZar)])
test3PZar<-as.data.frame(test2PZar)
frequbigPZar_lyrGenes=mergeByOverlaps(frequbigPZar_GRange,lyrGenes,type=c("any"))
frequbigPZar_lyrGenes_df=data.frame(as.data.frame(frequbigPZar_lyrGenes$frequbigPZar_GRange),as.data.frame(frequbigPZar_lyrGenes$lyrGenes),test3PZar$width)
colnames(frequbigPZar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigPZar_lyrGenes_df<-frequbigPZar_lyrGenes_df[(frequbigPZar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigPZar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigPZar_lyrGenes_df,"PZar_temp_WL300.xlsx",col.names=TRUE,row.names=FALSE)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigPZar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigPZar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigPZar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup=merge(frequbigPZar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PZar_WL300.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PZar_WL300.xlsx",overwrite=T)

Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigPZar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaPiekZapaarenosacnvregions_minRead6_WL300_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Piek arenosa#
#################################################

regionsZPar<-cnvr(resCNMOPSZPar)
singleZPar<-cnvs(resCNMOPSZPar)
testZPar<-findOverlaps(regionsZPar,singleZPar)
test2ZPar<- DataFrame(splitAsList(singleZPar$sampleName[subjectHits(testZPar)], queryHits(testZPar)),splitAsList(singleZPar$median[subjectHits(testZPar)], queryHits(testZPar)),splitAsList(singleZPar$mean[subjectHits(testZPar)], queryHits(testZPar)),splitAsList(singleZPar$CN[subjectHits(testZPar)], queryHits(testZPar)))
colnames(test2ZPar)<-c("sampleNames","mean","median","CN")
test7ZPar<-apply(as.data.frame(test2ZPar),1,dupl_cov)
mcols(regionsZPar)<-as.data.frame(as.matrix(test7ZPar))

dataZPar<-as.data.frame(regionsZPar)

adddataZPar<-do.call("rbind",apply(dataZPar,1,data_change))
dataZPar<-dataZPar[,-6]
frequZPar<-cbind(dataZPar,adddataZPar)

names(frequZPar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZPar<-frequZPar[count.fields(textConnection(frequZPar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZPar,"CnvregionswithmeanarenosaZapaPiekarenosa_minRead6_WL300_minL4_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZPar<-read.xlsx("CnvregionswithmeanarenosaZapaPiekarenosa_minRead6_WL300_minL4_min5samples.xlsx",1)
frequbigZPar[,2]<-as.numeric(as.character(frequbigZPar[,2]))
frequbigZPar[,3]<-as.numeric(as.character(frequbigZPar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZPar_GRange<-GRanges(seqnames=tolower(frequbigZPar$Scaffold),ranges=IRanges(start=frequbigZPar$start,end=frequbigZPar$end))
values(frequbigZPar_GRange)<-frequbigZPar[,6:11]

testZPar<-findOverlaps(frequbigZPar_GRange,lyrGenes)
test2ZPar<-pintersect(frequbigZPar_GRange[queryHits(testZPar)],lyrGenes[subjectHits(testZPar)])
test3ZPar<-as.data.frame(test2ZPar)
frequbigZPar_lyrGenes=mergeByOverlaps(frequbigZPar_GRange,lyrGenes,type=c("any"))
frequbigZPar_lyrGenes_df=data.frame(as.data.frame(frequbigZPar_lyrGenes$frequbigZPar_GRange),as.data.frame(frequbigZPar_lyrGenes$lyrGenes),test3ZPar$width)
colnames(frequbigZPar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZPar_lyrGenes_df<-frequbigZPar_lyrGenes_df[(frequbigZPar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZPar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZPar_lyrGenes_df,"ZPar_temp_WL300.xlsx",col.names=TRUE,row.names=FALSE)

frequbigZPar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZPar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZPar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup=merge(frequbigZPar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZPar_WL300.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZPar_WL300.xlsx",overwrite=T)

Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZPar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaZapaPiekarenosacnvregions_minRead6_WL300_minL4_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

PZar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PZar_WL300.xlsx",1)
ZPar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZPar_WL300.xlsx",1)

PZar_unique<-PZar[!duplicated(PZar$Lyr_Gene,PZar$Copy_start,PZar$Copy_end),]
ZPar_unique<-ZPar[!duplicated(ZPar$Lyr_Gene,ZPar$Copy_start,ZPar$Copy_end),]

PZar_unique$Lyr_Gene<-as.character(PZar_unique$Lyr_Gene)
ZPar_unique$Lyr_Gene<-as.character(ZPar_unique$Lyr_Gene)
PZar_unique$CN_class<-as.character(PZar_unique$CN_class)
ZPar_unique$CN_class<-as.character(ZPar_unique$CN_class)

PZaronly<-PZar_unique[!((PZar_unique$Lyr_Gene%in%ZPar_unique$Lyr_Gene)&(PZar_unique$CN_class==ZPar_unique$CN_class[match(PZar_unique$Lyr_Gene,ZPar_unique$Lyr_Gene)])),]
ZParonly<-ZPar_unique[!((ZPar_unique$Lyr_Gene%in%PZar_unique$Lyr_Gene)&(ZPar_unique$CN_class==PZar_unique$CN_class[match(ZPar_unique$Lyr_Gene,PZar_unique$Lyr_Gene)])),]

write.table(PZaronly,"PZar_cnvs_only_strict_WL300.table",row.names=F,sep="\t")
write.table(ZParonly,"ZPar_cnvs_only_strict_WL300.table",row.names=F,sep="\t")











