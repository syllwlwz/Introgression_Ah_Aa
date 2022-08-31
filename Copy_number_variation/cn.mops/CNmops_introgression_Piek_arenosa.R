#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesPiekarenosa <- scan("Piek_arenosa_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesPiekarenosa)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Piekarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesPiekarenosa,sampleNames=c("Piek_a10","Piek_a12","Piek_a13","Piek_a14","Piek_a15","Piek_a16","Piek_a17","Piek_a18","Piek_a19","Piek_a1","Piek_a20-1","Piek_a20-2","Piek_a22","Piek_a2","Piek_a4","Piek_a5","Piek_a7","Piek_a8","Piek_a9","Piek_h14_1","Piek_h6_2"),refSeqName=sort(scaffolds),WL=500) 
Kowaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesKowaarenosa,sampleNames=c("Kowa001a04","Kowa001a05","Kowa001a06","Kowa001a07","Kowa001a08","Kowa001a09","Kowa001a11","Kowa001a12"),refSeqName=sort(scaffolds),WL=500) 
PKar <- referencecn.mops(cases=Piekarenosa,controls=Kowaarenosa, minWidth = 2, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSPKar <- calcIntegerCopyNumbers(PKar)

KPar <- referencecn.mops(cases=Kowaarenosa,controls=Piekarenosa, minWidth = 2, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSKPar <- calcIntegerCopyNumbers(KPar)


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
regionsPKar<-cnvr(resCNMOPSPKar)
singlePKar<-cnvs(resCNMOPSPKar)
testPKar<-findOverlaps(regionsPKar,singlePKar)
test2PKar<- DataFrame(splitAsList(singlePKar$sampleName[subjectHits(testPKar)], queryHits(testPKar)),splitAsList(singlePKar$median[subjectHits(testPKar)], queryHits(testPKar)),splitAsList(singlePKar$mean[subjectHits(testPKar)], queryHits(testPKar)),splitAsList(singlePKar$CN[subjectHits(testPKar)], queryHits(testPKar)))
colnames(test2PKar)<-c("sampleNames","mean","median","CN")
test7PKar<-apply(as.data.frame(test2PKar),1,dupl_cov)
mcols(regionsPKar)<-as.data.frame(as.matrix(test7PKar))

dataPKar<-as.data.frame(regionsPKar)

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

adddataPKar<-do.call("rbind",apply(dataPKar,1,data_change))
dataPKar<-dataPKar[,-6]
frequPKar<-cbind(dataPKar,adddataPKar)

names(frequPKar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3PKar<-frequPKar[count.fields(textConnection(frequPKar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3PKar,"CnvregionswithmeanarenosaPiekKowaarenosa_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigPKar<-read.xlsx("CnvregionswithmeanarenosaPiekKowaarenosa_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigPKar[,2]<-as.numeric(as.character(frequbigPKar[,2]))
frequbigPKar[,3]<-as.numeric(as.character(frequbigPKar[,3]))

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
frequbigPKar_GRange<-GRanges(seqnames=tolower(frequbigPKar$Scaffold),ranges=IRanges(start=frequbigPKar$start,end=frequbigPKar$end))
values(frequbigPKar_GRange)<-frequbigPKar[,6:11]

testPKar<-findOverlaps(frequbigPKar_GRange,lyrGenes)
test2PKar<-pintersect(frequbigPKar_GRange[queryHits(testPKar)],lyrGenes[subjectHits(testPKar)])
test3PKar<-as.data.frame(test2PKar)
frequbigPKar_lyrGenes=mergeByOverlaps(frequbigPKar_GRange,lyrGenes,type=c("any"))
frequbigPKar_lyrGenes_df=data.frame(as.data.frame(frequbigPKar_lyrGenes$frequbigPKar_GRange),as.data.frame(frequbigPKar_lyrGenes$lyrGenes),test3PKar$width)
colnames(frequbigPKar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigPKar_lyrGenes_df<-frequbigPKar_lyrGenes_df[(frequbigPKar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigPKar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigPKar_lyrGenes_df,"PKar_temp.xlsx",col.names=TRUE,row.names=FALSE)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigPKar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigPKar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigPKar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup=merge(frequbigPKar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PKar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PKar.xlsx",overwrite=T)

Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigPKar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaPiekKowaarenosacnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Piek arenosa#
#################################################

regionsKPar<-cnvr(resCNMOPSKPar)
singleKPar<-cnvs(resCNMOPSKPar)
testKPar<-findOverlaps(regionsKPar,singleKPar)
test2KPar<- DataFrame(splitAsList(singleKPar$sampleName[subjectHits(testKPar)], queryHits(testKPar)),splitAsList(singleKPar$median[subjectHits(testKPar)], queryHits(testKPar)),splitAsList(singleKPar$mean[subjectHits(testKPar)], queryHits(testKPar)),splitAsList(singleKPar$CN[subjectHits(testKPar)], queryHits(testKPar)))
colnames(test2KPar)<-c("sampleNames","mean","median","CN")
test7KPar<-apply(as.data.frame(test2KPar),1,dupl_cov)
mcols(regionsKPar)<-as.data.frame(as.matrix(test7KPar))

dataKPar<-as.data.frame(regionsKPar)

adddataKPar<-do.call("rbind",apply(dataKPar,1,data_change))
dataKPar<-dataKPar[,-6]
frequKPar<-cbind(dataKPar,adddataKPar)

names(frequKPar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KPar<-frequKPar[count.fields(textConnection(frequKPar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KPar,"CnvregionswithmeanarenosaKowaPiekarenosa_minRead6_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKPar<-read.xlsx("CnvregionswithmeanarenosaKowaPiekarenosa_minRead6_WL500_minL2_min5samples.xlsx",1)
frequbigKPar[,2]<-as.numeric(as.character(frequbigKPar[,2]))
frequbigKPar[,3]<-as.numeric(as.character(frequbigKPar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKPar_GRange<-GRanges(seqnames=tolower(frequbigKPar$Scaffold),ranges=IRanges(start=frequbigKPar$start,end=frequbigKPar$end))
values(frequbigKPar_GRange)<-frequbigKPar[,6:11]

testKPar<-findOverlaps(frequbigKPar_GRange,lyrGenes)
test2KPar<-pintersect(frequbigKPar_GRange[queryHits(testKPar)],lyrGenes[subjectHits(testKPar)])
test3KPar<-as.data.frame(test2KPar)
frequbigKPar_lyrGenes=mergeByOverlaps(frequbigKPar_GRange,lyrGenes,type=c("any"))
frequbigKPar_lyrGenes_df=data.frame(as.data.frame(frequbigKPar_lyrGenes$frequbigKPar_GRange),as.data.frame(frequbigKPar_lyrGenes$lyrGenes),test3KPar$width)
colnames(frequbigKPar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKPar_lyrGenes_df<-frequbigKPar_lyrGenes_df[(frequbigKPar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKPar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKPar_lyrGenes_df,"KPar_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigKPar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKPar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKPar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup=merge(frequbigKPar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KPar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KPar.xlsx",overwrite=T)

Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKPar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaKowaPiekarenosacnvregions_minRead6_WL500_minL2_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

PKar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_PKar.xlsx",1)
KPar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KPar.xlsx",1)

PKar_unique<-PKar[!duplicated(PKar$Lyr_Gene,PKar$Copy_start,PKar$Copy_end),]
KPar_unique<-KPar[!duplicated(KPar$Lyr_Gene,KPar$Copy_start,KPar$Copy_end),]

PKar_unique$Lyr_Gene<-as.character(PKar_unique$Lyr_Gene)
KPar_unique$Lyr_Gene<-as.character(KPar_unique$Lyr_Gene)
PKar_unique$CN_class<-as.character(PKar_unique$CN_class)
KPar_unique$CN_class<-as.character(KPar_unique$CN_class)

PKaronly<-PKar_unique[!((PKar_unique$Lyr_Gene%in%KPar_unique$Lyr_Gene)&(PKar_unique$CN_class==KPar_unique$CN_class[match(PKar_unique$Lyr_Gene,KPar_unique$Lyr_Gene)])),]
KParonly<-KPar_unique[!((KPar_unique$Lyr_Gene%in%PKar_unique$Lyr_Gene)&(KPar_unique$CN_class==PKar_unique$CN_class[match(KPar_unique$Lyr_Gene,PKar_unique$Lyr_Gene)])),]

write.table(PKaronly,"PKar_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(KParonly,"KPar_cnvs_only_strict.table",row.names=F,sep="\t")











