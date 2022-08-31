#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesMiashalleri <- scan("Mias_halleri_full.list",character(),quote="")
BAMFilesLyrata <- scan("Lyrata_Vinod.list",character(),quote="")

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
test<-scanBamHeader(BAMFilesMiashalleri)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Miashalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesMiashalleri,sampleNames=c("Mias002h19","Mias002h20","Mias004h02","Mias004h06","Mias004h07","Mias004h08","Mias009h03","Mias009h05","Mias_a54","Mias_h09","Mias_h41","Mias_h43","Mias_h44","Mias_h45","Mias_h46"),refSeqName=sort(scaffolds),WL=350) 
Lyrata <- getReadCountsFromBAM(BAMFiles=BAMFilesLyrata,sampleNames=c("ICE1","Plech"),refSeqName=sort(scaffolds),WL=350)
MLha <- referencecn.mops(cases=Miashalleri,controls=Lyrata, minWidth = 6, minReadCount=6)
resCNMOPSMLha <- calcIntegerCopyNumbers(MLha)

LMha <- referencecn.mops(cases=Lyrata,controls=Miashalleri, minWidth = 6, minReadCount=6)
resCNMOPSLMha <- calcIntegerCopyNumbers(LMha)


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
regionsMLha<-cnvr(resCNMOPSMLha)
singleMLha<-cnvs(resCNMOPSMLha)
testMLha<-findOverlaps(regionsMLha,singleMLha)
test2MLha<- DataFrame(splitAsList(singleMLha$sampleName[subjectHits(testMLha)], queryHits(testMLha)),splitAsList(singleMLha$median[subjectHits(testMLha)], queryHits(testMLha)),splitAsList(singleMLha$mean[subjectHits(testMLha)], queryHits(testMLha)),splitAsList(singleMLha$CN[subjectHits(testMLha)], queryHits(testMLha)))
colnames(test2MLha)<-c("sampleNames","mean","median","CN")
test7MLha<-apply(as.data.frame(test2MLha),1,dupl_cov)
mcols(regionsMLha)<-as.data.frame(as.matrix(test7MLha))

dataMLha<-as.data.frame(regionsMLha)

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

adddataMLha<-do.call("rbind",apply(dataMLha,1,data_change))
dataMLha<-dataMLha[,-6]
frequMLha<-cbind(dataMLha,adddataMLha)

names(frequMLha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3MLha<-frequMLha[count.fields(textConnection(frequMLha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3MLha,"CnvregionswithmeanhalleriMiasLyratahalleri_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigMLha<-read.xlsx("CnvregionswithmeanhalleriMiasLyratahalleri_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigMLha[,2]<-as.numeric(as.character(frequbigMLha[,2]))
frequbigMLha[,3]<-as.numeric(as.character(frequbigMLha[,3]))

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
frequbigMLha_GRange<-GRanges(seqnames=tolower(frequbigMLha$Scaffold),ranges=IRanges(start=frequbigMLha$start,end=frequbigMLha$end))
values(frequbigMLha_GRange)<-frequbigMLha[,6:11]

testMLha<-findOverlaps(frequbigMLha_GRange,lyrGenes)
test2MLha<-pintersect(frequbigMLha_GRange[queryHits(testMLha)],lyrGenes[subjectHits(testMLha)])
test3MLha<-as.data.frame(test2MLha)
frequbigMLha_lyrGenes=mergeByOverlaps(frequbigMLha_GRange,lyrGenes,type=c("any"))
frequbigMLha_lyrGenes_df=data.frame(as.data.frame(frequbigMLha_lyrGenes$frequbigMLha_GRange),as.data.frame(frequbigMLha_lyrGenes$lyrGenes),test3MLha$width)
colnames(frequbigMLha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigMLha_lyrGenes_df<-frequbigMLha_lyrGenes_df[(frequbigMLha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigMLha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigMLha_lyrGenes_df,"MLha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigMLha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigMLha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigMLha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup=merge(frequbigMLha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MLha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MLha.xlsx",overwrite=T)

Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigMLha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriMiasLyratahallericnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Lyrata to Mias halleri#
#################################################

regionsLMha<-cnvr(resCNMOPSLMha)
singleLMha<-cnvs(resCNMOPSLMha)
testLMha<-findOverlaps(regionsLMha,singleLMha)
test2LMha<- DataFrame(splitAsList(singleLMha$sampleName[subjectHits(testLMha)], queryHits(testLMha)),splitAsList(singleLMha$median[subjectHits(testLMha)], queryHits(testLMha)),splitAsList(singleLMha$mean[subjectHits(testLMha)], queryHits(testLMha)),splitAsList(singleLMha$CN[subjectHits(testLMha)], queryHits(testLMha)))
colnames(test2LMha)<-c("sampleNames","mean","median","CN")
test7LMha<-apply(as.data.frame(test2LMha),1,dupl_cov)
mcols(regionsLMha)<-as.data.frame(as.matrix(test7LMha))

dataLMha<-as.data.frame(regionsLMha)

adddataLMha<-do.call("rbind",apply(dataLMha,1,data_change))
dataLMha<-dataLMha[,-6]
frequLMha<-cbind(dataLMha,adddataLMha)

names(frequLMha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3LMha<-frequLMha[count.fields(textConnection(frequLMha$sampleNames),sep=",")>=2,]

require(openxlsx)
write.xlsx(frequ3LMha,"CnvregionswithmeanhalleriLyrataMiashalleri_minRead6_WL350_minL6_min2samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigLMha<-read.xlsx("CnvregionswithmeanhalleriLyrataMiashalleri_minRead6_WL350_minL6_min2samples.xlsx",1)
frequbigLMha[,2]<-as.numeric(as.character(frequbigLMha[,2]))
frequbigLMha[,3]<-as.numeric(as.character(frequbigLMha[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigLMha_GRange<-GRanges(seqnames=tolower(frequbigLMha$Scaffold),ranges=IRanges(start=frequbigLMha$start,end=frequbigLMha$end))
values(frequbigLMha_GRange)<-frequbigLMha[,6:11]

testLMha<-findOverlaps(frequbigLMha_GRange,lyrGenes)
test2LMha<-pintersect(frequbigLMha_GRange[queryHits(testLMha)],lyrGenes[subjectHits(testLMha)])
test3LMha<-as.data.frame(test2LMha)
frequbigLMha_lyrGenes=mergeByOverlaps(frequbigLMha_GRange,lyrGenes,type=c("any"))
frequbigLMha_lyrGenes_df=data.frame(as.data.frame(frequbigLMha_lyrGenes$frequbigLMha_GRange),as.data.frame(frequbigLMha_lyrGenes$lyrGenes),test3LMha$width)
colnames(frequbigLMha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigLMha_lyrGenes_df<-frequbigLMha_lyrGenes_df[(frequbigLMha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigLMha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigLMha_lyrGenes_df,"LMha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigLMha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigLMha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigLMha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup=merge(frequbigLMha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LMha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LMha.xlsx",overwrite=T)

Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigLMha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriLyrataMiashallericnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast2samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

MLha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_MLha.xlsx",1)
LMha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_LMha.xlsx",1)

MLha_unique<-MLha[!duplicated(MLha$Lyr_Gene,MLha$Copy_start,MLha$Copy_end),]
LMha_unique<-LMha[!duplicated(LMha$Lyr_Gene,LMha$Copy_start,LMha$Copy_end),]

MLha_unique$Lyr_Gene<-as.character(MLha_unique$Lyr_Gene)
LMha_unique$Lyr_Gene<-as.character(LMha_unique$Lyr_Gene)
MLha_unique$CN_class<-as.character(MLha_unique$CN_class)
LMha_unique$CN_class<-as.character(LMha_unique$CN_class)

MLhaonly<-MLha_unique[!((MLha_unique$Lyr_Gene%in%LMha_unique$Lyr_Gene)&(MLha_unique$CN_class==LMha_unique$CN_class[match(MLha_unique$Lyr_Gene,LMha_unique$Lyr_Gene)])),]
LMhaonly<-LMha_unique[!((LMha_unique$Lyr_Gene%in%MLha_unique$Lyr_Gene)&(LMha_unique$CN_class==MLha_unique$CN_class[match(LMha_unique$Lyr_Gene,MLha_unique$Lyr_Gene)])),]

write.table(MLhaonly,"MLha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(LMhaonly,"LMha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











