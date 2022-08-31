#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesKatoarenosa <- scan("Kato_arenosa_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesKatoarenosa)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Katoarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesKatoarenosa,sampleNames=c("Kato_a21","Kato_a22","Kato_a23","Kato_a24","Kato_a27","Kato_a28","Kato_a29","Kato_a30","Kato_a33","Kato_a35","Kato_h21","Kato_h22","Kato_h23","Kato_h24","Kato_h26","Kato_h27","Kato_h33","Kato_h35"),refSeqName=sort(scaffolds),WL=350) 
Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapaarenosa,sampleNames=c("Zapa002a01","Zapa002a02","Zapa002a03","Zapa002a04","Zapa002a05","Zapa002a07","Zapa004a11","Zapa008a19","Zapa_11a03x12a01_l","Zapa_a09xa03_b"),refSeqName=sort(scaffolds),WL=350)
KZar <- referencecn.mops(cases=Katoarenosa,controls=Zapaarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSKZar <- calcIntegerCopyNumbers(KZar)

ZKar <- referencecn.mops(cases=Zapaarenosa,controls=Katoarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSZKar <- calcIntegerCopyNumbers(ZKar)


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
regionsKZar<-cnvr(resCNMOPSKZar)
singleKZar<-cnvs(resCNMOPSKZar)
testKZar<-findOverlaps(regionsKZar,singleKZar)
test2KZar<- DataFrame(splitAsList(singleKZar$sampleName[subjectHits(testKZar)], queryHits(testKZar)),splitAsList(singleKZar$median[subjectHits(testKZar)], queryHits(testKZar)),splitAsList(singleKZar$mean[subjectHits(testKZar)], queryHits(testKZar)),splitAsList(singleKZar$CN[subjectHits(testKZar)], queryHits(testKZar)))
colnames(test2KZar)<-c("sampleNames","mean","median","CN")
test7KZar<-apply(as.data.frame(test2KZar),1,dupl_cov)
mcols(regionsKZar)<-as.data.frame(as.matrix(test7KZar))

dataKZar<-as.data.frame(regionsKZar)

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

adddataKZar<-do.call("rbind",apply(dataKZar,1,data_change))
dataKZar<-dataKZar[,-6]
frequKZar<-cbind(dataKZar,adddataKZar)

names(frequKZar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KZar<-frequKZar[count.fields(textConnection(frequKZar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KZar,"CnvregionswithmeanarenosaKatoZapaarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigKZar<-read.xlsx("CnvregionswithmeanarenosaKatoZapaarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigKZar[,2]<-as.numeric(as.character(frequbigKZar[,2]))
frequbigKZar[,3]<-as.numeric(as.character(frequbigKZar[,3]))

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
frequbigKZar_GRange<-GRanges(seqnames=tolower(frequbigKZar$Scaffold),ranges=IRanges(start=frequbigKZar$start,end=frequbigKZar$end))
values(frequbigKZar_GRange)<-frequbigKZar[,6:11]

testKZar<-findOverlaps(frequbigKZar_GRange,lyrGenes)
test2KZar<-pintersect(frequbigKZar_GRange[queryHits(testKZar)],lyrGenes[subjectHits(testKZar)])
test3KZar<-as.data.frame(test2KZar)
frequbigKZar_lyrGenes=mergeByOverlaps(frequbigKZar_GRange,lyrGenes,type=c("any"))
frequbigKZar_lyrGenes_df=data.frame(as.data.frame(frequbigKZar_lyrGenes$frequbigKZar_GRange),as.data.frame(frequbigKZar_lyrGenes$lyrGenes),test3KZar$width)
colnames(frequbigKZar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKZar_lyrGenes_df<-frequbigKZar_lyrGenes_df[(frequbigKZar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKZar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKZar_lyrGenes_df,"KaZar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigKZar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKZar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKZar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup=merge(frequbigKZar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KaZar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KaZar.xlsx",overwrite=T)

Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigKZar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaKatoZapaarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to Kato arenosa#
#################################################

regionsZKar<-cnvr(resCNMOPSZKar)
singleZKar<-cnvs(resCNMOPSZKar)
testZKar<-findOverlaps(regionsZKar,singleZKar)
test2ZKar<- DataFrame(splitAsList(singleZKar$sampleName[subjectHits(testZKar)], queryHits(testZKar)),splitAsList(singleZKar$median[subjectHits(testZKar)], queryHits(testZKar)),splitAsList(singleZKar$mean[subjectHits(testZKar)], queryHits(testZKar)),splitAsList(singleZKar$CN[subjectHits(testZKar)], queryHits(testZKar)))
colnames(test2ZKar)<-c("sampleNames","mean","median","CN")
test7ZKar<-apply(as.data.frame(test2ZKar),1,dupl_cov)
mcols(regionsZKar)<-as.data.frame(as.matrix(test7ZKar))

dataZKar<-as.data.frame(regionsZKar)

adddataZKar<-do.call("rbind",apply(dataZKar,1,data_change))
dataZKar<-dataZKar[,-6]
frequZKar<-cbind(dataZKar,adddataZKar)

names(frequZKar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZKar<-frequZKar[count.fields(textConnection(frequZKar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZKar,"CnvregionswithmeanarenosaZapaKatoarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZKar<-read.xlsx("CnvregionswithmeanarenosaZapaKatoarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigZKar[,2]<-as.numeric(as.character(frequbigZKar[,2]))
frequbigZKar[,3]<-as.numeric(as.character(frequbigZKar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZKar_GRange<-GRanges(seqnames=tolower(frequbigZKar$Scaffold),ranges=IRanges(start=frequbigZKar$start,end=frequbigZKar$end))
values(frequbigZKar_GRange)<-frequbigZKar[,6:11]

testZKar<-findOverlaps(frequbigZKar_GRange,lyrGenes)
test2ZKar<-pintersect(frequbigZKar_GRange[queryHits(testZKar)],lyrGenes[subjectHits(testZKar)])
test3ZKar<-as.data.frame(test2ZKar)
frequbigZKar_lyrGenes=mergeByOverlaps(frequbigZKar_GRange,lyrGenes,type=c("any"))
frequbigZKar_lyrGenes_df=data.frame(as.data.frame(frequbigZKar_lyrGenes$frequbigZKar_GRange),as.data.frame(frequbigZKar_lyrGenes$lyrGenes),test3ZKar$width)
colnames(frequbigZKar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZKar_lyrGenes_df<-frequbigZKar_lyrGenes_df[(frequbigZKar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZKar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZKar_lyrGenes_df,"ZKaar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZKar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZKar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZKar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup=merge(frequbigZKar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZKaar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZKaar.xlsx",overwrite=T)

Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZKar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaZapaKatoarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="GenesarenosaWulmoutgroupcnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

KZar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KaZar.xlsx",1)
ZKar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZKaar.xlsx",1)

KZar_unique<-KZar[!duplicated(KZar$Lyr_Gene,KZar$Copy_start,KZar$Copy_end),]
ZKar_unique<-ZKar[!duplicated(ZKar$Lyr_Gene,ZKar$Copy_start,ZKar$Copy_end),]

KZar_unique$Lyr_Gene<-as.character(KZar_unique$Lyr_Gene)
ZKar_unique$Lyr_Gene<-as.character(ZKar_unique$Lyr_Gene)
KZar_unique$CN_class<-as.character(KZar_unique$CN_class)
ZKar_unique$CN_class<-as.character(ZKar_unique$CN_class)

KZaronly<-KZar_unique[!((KZar_unique$Lyr_Gene%in%ZKar_unique$Lyr_Gene)&(KZar_unique$CN_class==ZKar_unique$CN_class[match(KZar_unique$Lyr_Gene,ZKar_unique$Lyr_Gene)])),]
ZKaronly<-ZKar_unique[!((ZKar_unique$Lyr_Gene%in%KZar_unique$Lyr_Gene)&(ZKar_unique$CN_class==KZar_unique$CN_class[match(ZKar_unique$Lyr_Gene,KZar_unique$Lyr_Gene)])),]

write.table(KZaronly,"KaZar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(ZKaronly,"ZKaar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











