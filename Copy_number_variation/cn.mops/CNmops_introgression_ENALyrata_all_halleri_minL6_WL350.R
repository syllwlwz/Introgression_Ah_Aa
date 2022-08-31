#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesENALyrata <- scan("ENALyrata.list",character(),quote="")
BAMFileshalleri <- scan("All_halleri_full.list",character(),quote="")

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
halleri <- getReadCountsFromBAM(BAMFiles=BAMFileshalleri,sampleNames=c("Kato_h09","Buko_h1_7","Buko_h04","Buko_h12","Buko_h21","Buko_h22","Buko_h23","Buko_h24","Buko_h25","Buko_h26","Buko_h27","Buko_h28","Buko_h29","Buko_h30","Buko_h31","Buko_h32","Buko_h33",
"Buko_h35","Buko_h36","Mias_a54","Mias_h09","Mias_h41","Mias_h43","Mias_h44","Mias_h45","Mias_h46","Mias002h19","Mias002h20","Mias004h02","Mias004h06","Mias004h07","Mias004h08","Mias009h03","Mias009h05","Piek_h01","Piek_h03","Piek_h04","Piek_h05comb","
Piek_h06comb","Piek_h07","Piek_h10","Piek_h11","Piek_h12","Piek_h13","Piek_h14_2","Piek_h15","Piek_h2","Kowa002h13","Kowa002h14","Kowa002h17","Kowa002h18","Kowa002h21","Kowa_h02comb","Kowa_h04comb","Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43","Zako_h01comb","
Zako_h03","Zako_h09comb","Zako_h10","Zako_h12","Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zapa_h01comb","Zapa_h07","Zapa_h11comb","Zapa_h12"),refSeqName=sort(scaffolds),WL=350) 
ENALha <- referencecn.mops(cases=ENALyrata,controls=halleri, minWidth = 6, minReadCount=6)
resCNMOPSENALha <- calcIntegerCopyNumbers(ENALha)

#ENALyrataarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesENALyrataarenosa,sampleNames=c("ENALyrata_h04","ENALyrata_h12","ENALyrata_h1_7","ENALyrata_h21","ENALyrata_h22","ENALyrata_h23","ENALyrata_h24","ENALyrata_h25","ENALyrata_h26","ENALyrata_h27","ENALyrata_h28","ENALyrata_h29","ENALyrata_h30","ENALyrata_h31","ENALyrata_h32","ENALyrata_h33","ENALyrata_h35","ENALyrata_h36"),refSeqName=sort(scaffolds),WL=500) 
#Zapaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesZapahalleri,sampleNames=c("Zapa002h13","Zapa002h14","Zapa002h17","Zapa002h18","Zapa002h21","Zapa_h02comb","Zapa_h04comb","Zapa_h07","Zapa_h41","Zapa_h42","Zapa_h43"),refSeqName=sort(scaffolds),WL=500) 
#BKa <- referencecn.mops(cases=ENALyrataarenosa,controls=Zapaarenosa, minWidth = 2, minReadCount=6)
#resCNMOPSENALha <- calcIntegerCopyNumbers(BKa)


ZENALha <- referencecn.mops(cases=halleri,controls=ENALyrata, minWidth = 6, minReadCount=6)
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
regionsENALha<-cnvr(resCNMOPSENALha)
singleENALha<-cnvs(resCNMOPSENALha)
testENALha<-findOverlaps(regionsENALha,singleENALha)
test2ENALha<- DataFrame(splitAsList(singleENALha$sampleName[subjectHits(testENALha)], queryHits(testENALha)),splitAsList(singleENALha$median[subjectHits(testENALha)], queryHits(testENALha)),splitAsList(singleENALha$mean[subjectHits(testENALha)], queryHits(testENALha)),splitAsList(singleENALha$CN[subjectHits(testENALha)], queryHits(testENALha)))
colnames(test2ENALha)<-c("sampleNames","mean","median","CN")
test7ENALha<-apply(as.data.frame(test2ENALha),1,dupl_cov)
mcols(regionsENALha)<-as.data.frame(as.matrix(test7ENALha))

dataENALha<-as.data.frame(regionsENALha)

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

adddataENALha<-do.call("rbind",apply(dataENALha,1,data_change))
dataENALha<-dataENALha[,-6]
frequENALha<-cbind(dataENALha,adddataENALha)

names(frequENALha)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3ENALha<-frequENALha[count.fields(textConnection(frequENALha$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3ENALha,"CnvregionswithmeanhalleriENALyrataZapahalleri_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigENALha<-read.xlsx("CnvregionswithmeanhalleriENALyrataZapahalleri_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigENALha[,2]<-as.numeric(as.character(frequbigENALha[,2]))
frequbigENALha[,3]<-as.numeric(as.character(frequbigENALha[,3]))

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
frequbigENALha_GRange<-GRanges(seqnames=tolower(frequbigENALha$Scaffold),ranges=IRanges(start=frequbigENALha$start,end=frequbigENALha$end))
values(frequbigENALha_GRange)<-frequbigENALha[,6:11]

testENALha<-findOverlaps(frequbigENALha_GRange,lyrGenes)
test2ENALha<-pintersect(frequbigENALha_GRange[queryHits(testENALha)],lyrGenes[subjectHits(testENALha)])
test3ENALha<-as.data.frame(test2ENALha)
frequbigENALha_lyrGenes=mergeByOverlaps(frequbigENALha_GRange,lyrGenes,type=c("any"))
frequbigENALha_lyrGenes_df=data.frame(as.data.frame(frequbigENALha_lyrGenes$frequbigENALha_GRange),as.data.frame(frequbigENALha_lyrGenes$lyrGenes),test3ENALha$width)
colnames(frequbigENALha_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigENALha_lyrGenes_df<-frequbigENALha_lyrGenes_df[(frequbigENALha_lyrGenes_df$Widthofoverlap>=(0.9*frequbigENALha_lyrGenes_df$gene_size)),]
#write.xlsx(frequbigENALha_lyrGenes_df,"ENALha_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigENALha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Alyr_ID",all.x=TRUE)
colnames(frequbigENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup=merge(frequbigENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ENALha.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ENALha.xlsx",overwrite=T)

Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigENALha_lyrGenes_df_desc_w_orthogroup_highsupport,"GeneshalleriENALyrataallhallericnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

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
write.xlsx(frequ3ZENALha,"CnvregionswithmeanhalleriENALyrata_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigZENALha<-read.xlsx("CnvregionswithmeanhalleriENALyrata_minRead6_WL350_minL6_min5samples.xlsx",1)
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

frequbigZENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigZENALha_lyrGenes_df,Thal_description, by.x="Gene", by.y="Alyr_ID",all.x=TRUE)
colnames(frequbigZENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup=merge(frequbigZENALha_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_haENAL.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_haENAL.xlsx",overwrite=T)

Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigZENALha_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesallhalleriENALyratacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Geneshallericnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

ENALha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_ENALha.xlsx",1)
ZENALha<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_haENAL.xlsx",1)

ENALha_unique<-ENALha[!duplicated(ENALha$Lyr_Gene,ENALha$Copy_start,ENALha$Copy_end),]
ZENALha_unique<-ZENALha[!duplicated(ZENALha$Lyr_Gene,ZENALha$Copy_start,ZENALha$Copy_end),]

ENALha_unique$Lyr_Gene<-as.character(ENALha_unique$Lyr_Gene)
ZENALha_unique$Lyr_Gene<-as.character(ZENALha_unique$Lyr_Gene)
ENALha_unique$CN_class<-as.character(ENALha_unique$CN_class)
ZENALha_unique$CN_class<-as.character(ZENALha_unique$CN_class)

ENALhaonly<-ENALha_unique[!((ENALha_unique$Lyr_Gene%in%ZENALha_unique$Lyr_Gene)&(ENALha_unique$CN_class==ZENALha_unique$CN_class[match(ENALha_unique$Lyr_Gene,ZENALha_unique$Lyr_Gene)])),]
ZENALhaonly<-ZENALha_unique[!((ZENALha_unique$Lyr_Gene%in%ENALha_unique$Lyr_Gene)&(ZENALha_unique$CN_class==ENALha_unique$CN_class[match(ZENALha_unique$Lyr_Gene,ENALha_unique$Lyr_Gene)])),]

write.table(ENALhaonly,"ENALha_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(ZENALhaonly,"haENAL_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











