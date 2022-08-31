#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesENALyrata <- scan("ENALyrata.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesENALyrata)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

ENALyrata <- getReadCountsFromBAM(BAMFiles=BAMFilesENALyrata,sampleNames=c("A_ly_1","A_ly_2","A_ly_3","A_ly_4","A_ly_5","A_ly_6","A_ly_7"),refSeqName=sort(scaffolds),WL=350) 
Zapahalleri <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zako_h01comb","Zako_h03","Zako_h09comb","Zako_h10","Zako_h12","Zapa_h01comb","Zapa_h07","Zapa_h11comb","Zapa_h12"),refSeqName=sort(scaffolds),WL=350) 
ENALZha <- referencecn.mops(cases=ENALyrata,controls=Zapahalleri, minWidth = 6, minReadCount=6)
resCNMOPSENALZha <- calcIntegerCopyNumbers(ENALZha)

#ENALyrataarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesENALyrataarenosa,sampleNames=c("ENALyrata_h04","ENALyrata_h12","ENALyrata_h1_7","ENALyrata_h21","ENALyrata_h22","ENALyrata_h23","ENALyrata_h24","ENALyrata_h25","ENALyrata_h26","ENALyrata_h27","ENALyrata_h28","ENALyrata_h29","ENALyrata_h30","ENALyrata_h31","ENALyrata_h32","ENALyrata_h33","ENALyrata_h35","ENALyrata_h36"),refSeqName=sort(scaffolds),WL=500) 
#Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zapa002h13","Zapa002h14","Zapa002h17","Zapa002h18","Zapa002h21","Zapa_h02comb","Zapa_h04comb","Zapa_h07","Zapa_h41","Zapa_h42","Zapa_h43"),refSeqName=sort(scaffolds),WL=500) 
#BKa <- referencecn.mops(cases=ENALyrataarenosa,controls=Zapaarenosa, minWidth = 2, minReadCount=6)
#resCNMOPSENALZha <- calcIntegerCopyNumbers(BKa)


ZENALha <- referencecn.mops(cases=Zapahalleri,controls=ENALyrata, minWidth = 6, minReadCount=6)
resCNMOPSZENALha <- calcIntegerCopyNumbers(ZENALha)


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
regionsENALZha<-cnvr(resCNMOPSENALZha)
singleENALZha<-cnvs(resCNMOPSENALZha)
testENALZha<-findOverlaps(regionsENALZha,singleENALZha)
test2ENALZha<- DataFrame(splitAsList(singleENALZha$sampleName[subjectHits(testENALZha)], queryHits(testENALZha)),splitAsList(singleENALZha$median[subjectHits(testENALZha)], queryHits(testENALZha)),splitAsList(singleENALZha$mean[subjectHits(testENALZha)], queryHits(testENALZha)),splitAsList(singleENALZha$CN[subjectHits(testENALZha)], queryHits(testENALZha)))
colnames(test2ENALZha)<-c("sampleNames","mean","median","CN")
test7ENALZha<-apply(as.data.frame(test2ENALZha),1,dupl_cov)
mcols(regionsENALZha)<-as.data.frame(as.matrix(test7ENALZha))

dataENALZha<-as.data.frame(regionsENALZha)

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

adddataENALZha<-do.call("rbind",apply(dataENALZha,1,data_change))
dataENALZha<-dataENALZha[,-6]
frequENALZha<-cbind(dataENALZha,adddataENALZha)

names(frequENALZha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ENALZha<-frequENALZha[count.fields(textConnection(frequENALZha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ENALZha,"CnvregionswithmeanhalleriENALyrataZapahalleri_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigENALZha<-read.xlsx("CnvregionswithmeanhalleriENALyrataZapahalleri_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigENALZha[,2]<-as.numeric(as.character(frequbigENALZha[,2]))
frequbigENALZha[,3]<-as.numeric(as.character(frequbigENALZha[,3]))

###################################################
#find overlaps between windows and A. ENALyrata genes#
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
frequbigENALZha_GRange<-GRanges(seqnames=tolower(frequbigENALZha$Scaffold),ranges=IRanges(start=frequbigENALZha$start,end=frequbigENALZha$end))
values(frequbigENALZha_GRange)<-frequbigENALZha[,6:11]

testENALZha<-findOverlaps(frequbigENALZha_GRange,lyrGenes)
test2ENALZha<-pintersect(frequbigENALZha_GRange[queryHits(testENALZha)],lyrGenes[subjectHits(testENALZha)])
test3ENALZha<-as.data.frame(test2ENALZha)
frequbigENALZha_lyrGenes=mergeByOverlaps(frequbigENALZha_GRange,lyrGenes,type=c("any"))
frequbigENALZha_lyrGenes_df=data.frame(as.data.frame(frequbigENALZha_lyrGenes$frequbigENALZha_GRange),as.data.frame(frequbigENALZha_lyrGenes$lyrGenes),test3ENALZha$width)
colnames(frequbigENALZha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigENALZha_lyrGenes_df<-frequbigENALZha_lyrGenes_df[(frequbigENALZha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigENALZha_lyrGenes_df$gene_size)),]
#write.xlsx(frequbigENALZha_lyrGenes_df,"ENALZha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigENALZha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigENALZha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigENALZha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup=merge(frequbigENALZha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ENALZha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ENALZha.xlsx",overwrite=T)

Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigENALZha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriENALyrataZapahallericnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Zapa to ENALyrata halleri#
#################################################

regionsZENALha<-cnvr(resCNMOPSZENALha)
singleZENALha<-cnvs(resCNMOPSZENALha)
testZENALha<-findOverlaps(regionsZENALha,singleZENALha)
test2ZENALha<- DataFrame(splitAsList(singleZENALha$sampleName[subjectHits(testZENALha)], queryHits(testZENALha)),splitAsList(singleZENALha$median[subjectHits(testZENALha)], queryHits(testZENALha)),splitAsList(singleZENALha$mean[subjectHits(testZENALha)], queryHits(testZENALha)),splitAsList(singleZENALha$CN[subjectHits(testZENALha)], queryHits(testZENALha)))
colnames(test2ZENALha)<-c("sampleNames","mean","median","CN")
test7ZENALha<-apply(as.data.frame(test2ZENALha),1,dupl_cov)
mcols(regionsZENALha)<-as.data.frame(as.matrix(test7ZENALha))

dataZENALha<-as.data.frame(regionsZENALha)

adddataZENALha<-do.call("rbind",apply(dataZENALha,1,data_change))
dataZENALha<-dataZENALha[,-6]
frequZENALha<-cbind(dataZENALha,adddataZENALha)

names(frequZENALha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ZENALha<-frequZENALha[count.fields(textConnection(frequZENALha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ZENALha,"CnvregionswithmeanhalleriZapaENALyrata_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZENALha<-read.xlsx("CnvregionswithmeanhalleriZapaENALyrata_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigZENALha[,2]<-as.numeric(as.character(frequbigZENALha[,2]))
frequbigZENALha[,3]<-as.numeric(as.character(frequbigZENALha[,3]))

####################################################
#find overlaps between windows and A. ENALyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigZENALha_GRange<-GRanges(seqnames=tolower(frequbigZENALha$Scaffold),ranges=IRanges(start=frequbigZENALha$start,end=frequbigZENALha$end))
values(frequbigZENALha_GRange)<-frequbigZENALha[,6:11]

testZENALha<-findOverlaps(frequbigZENALha_GRange,lyrGenes)
test2ZENALha<-pintersect(frequbigZENALha_GRange[queryHits(testZENALha)],lyrGenes[subjectHits(testZENALha)])
test3ZENALha<-as.data.frame(test2ZENALha)
frequbigZENALha_lyrGenes=mergeByOverlaps(frequbigZENALha_GRange,lyrGenes,type=c("any"))
frequbigZENALha_lyrGenes_df=data.frame(as.data.frame(frequbigZENALha_lyrGenes$frequbigZENALha_GRange),as.data.frame(frequbigZENALha_lyrGenes$lyrGenes),test3ZENALha$width)
colnames(frequbigZENALha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigZENALha_lyrGenes_df<-frequbigZENALha_lyrGenes_df[(frequbigZENALha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigZENALha_lyrGenes_df$gene_size)),]
write.xlsx(frequbigZENALha_lyrGenes_df,"ZENALha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigZENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZENALha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigZENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup=merge(frequbigZENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZENALha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZENALha.xlsx",overwrite=T)

Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriZapaENALyratacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

ENALZha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ENALZha.xlsx",1)
ZENALha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ZENALha.xlsx",1)

ENALZha_unique<-ENALZha[!duplicated(ENALZha$Lyr_Gene,ENALZha$Copy_start,ENALZha$Copy_end),]
ZENALha_unique<-ZENALha[!duplicated(ZENALha$Lyr_Gene,ZENALha$Copy_start,ZENALha$Copy_end),]

ENALZha_unique$Lyr_Gene<-as.character(ENALZha_unique$Lyr_Gene)
ZENALha_unique$Lyr_Gene<-as.character(ZENALha_unique$Lyr_Gene)
ENALZha_unique$CN_class<-as.character(ENALZha_unique$CN_class)
ZENALha_unique$CN_class<-as.character(ZENALha_unique$CN_class)

ENALZhaonly<-ENALZha_unique[!((ENALZha_unique$Lyr_Gene%in%ZENALha_unique$Lyr_Gene)&(ENALZha_unique$CN_class==ZENALha_unique$CN_class[match(ENALZha_unique$Lyr_Gene,ZENALha_unique$Lyr_Gene)])),]
ZENALhaonly<-ZENALha_unique[!((ZENALha_unique$Lyr_Gene%in%ENALZha_unique$Lyr_Gene)&(ZENALha_unique$CN_class==ENALZha_unique$CN_class[match(ZENALha_unique$Lyr_Gene,ENALZha_unique$Lyr_Gene)])),]

write.table(ENALZhaonly,"ENALZha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(ZENALhaonly,"ZENALha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











