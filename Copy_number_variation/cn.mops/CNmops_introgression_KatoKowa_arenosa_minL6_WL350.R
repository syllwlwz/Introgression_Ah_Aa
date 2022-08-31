#Cnmops
options(java.parameters = "-Xmx12000m")

require(cn.mops)

BAMFilesKatoarenosa <- scan("Kato_arenosa_full.list",character(),quote="")
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
test<-scanBamHeader(BAMFilesKatoarenosa)
scaffolds<-names(test[[1]]$targets)[1:9]

require(DNAcopy)

Katoarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesKatoarenosa,sampleNames=c("Kato_a21","Kato_a22","Kato_a23","Kato_a24","Kato_a27","Kato_a28","Kato_a29","Kato_a30","Kato_a33","Kato_a35","Kato_h21","Kato_h22","Kato_h23","Kato_h24","Kato_h26","Kato_h27","Kato_h33","Kato_h35"),refSeqName=sort(scaffolds),WL=350) 
Kowaarenosa <- getReadCountsFromBAM(BAMFiles=BAMFilesKowaarenosa,sampleNames=c("Kowa001a04","Kowa001a05","Kowa001a06","Kowa001a07","Kowa001a08","Kowa001a09","Kowa001a11","Kowa001a12"),refSeqName=sort(scaffolds),WL=350) 
KaKoar <- referencecn.mops(cases=Katoarenosa,controls=Kowaarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSKaKoar <- calcIntegerCopyNumbers(KaKoar)

KKar <- referencecn.mops(cases=Kowaarenosa,controls=Katoarenosa, minWidth = 6, minReadCount=6,I=c(0.0125,0.025,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,4),classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN16"))
resCNMOPSKKar <- calcIntegerCopyNumbers(KKar)


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
regionsKaKoar<-cnvr(resCNMOPSKaKoar)
singleKaKoar<-cnvs(resCNMOPSKaKoar)
testKaKoar<-findOverlaps(regionsKaKoar,singleKaKoar)
test2KaKoar<- DataFrame(splitAsList(singleKaKoar$sampleName[subjectHits(testKaKoar)], queryHits(testKaKoar)),splitAsList(singleKaKoar$median[subjectHits(testKaKoar)], queryHits(testKaKoar)),splitAsList(singleKaKoar$mean[subjectHits(testKaKoar)], queryHits(testKaKoar)),splitAsList(singleKaKoar$CN[subjectHits(testKaKoar)], queryHits(testKaKoar)))
colnames(test2KaKoar)<-c("sampleNames","mean","median","CN")
test7KaKoar<-apply(as.data.frame(test2KaKoar),1,dupl_cov)
mcols(regionsKaKoar)<-as.data.frame(as.matrix(test7KaKoar))

dataKaKoar<-as.data.frame(regionsKaKoar)

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

adddataKaKoar<-do.call("rbind",apply(dataKaKoar,1,data_change))
dataKaKoar<-dataKaKoar[,-6]
frequKaKoar<-cbind(dataKaKoar,adddataKaKoar)

names(frequKaKoar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KaKoar<-frequKaKoar[count.fields(textConnection(frequKaKoar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KaKoar,"CnvregionswithmeanarenosaKatoKowaarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigKaKoar<-read.xlsx("CnvregionswithmeanarenosaKatoKowaarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigKaKoar[,2]<-as.numeric(as.character(frequbigKaKoar[,2]))
frequbigKaKoar[,3]<-as.numeric(as.character(frequbigKaKoar[,3]))

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
frequbigKaKoar_GRange<-GRanges(seqnames=tolower(frequbigKaKoar$Scaffold),ranges=IRanges(start=frequbigKaKoar$start,end=frequbigKaKoar$end))
values(frequbigKaKoar_GRange)<-frequbigKaKoar[,6:11]

testKaKoar<-findOverlaps(frequbigKaKoar_GRange,lyrGenes)
test2KaKoar<-pintersect(frequbigKaKoar_GRange[queryHits(testKaKoar)],lyrGenes[subjectHits(testKaKoar)])
test3KaKoar<-as.data.frame(test2KaKoar)
frequbigKaKoar_lyrGenes=mergeByOverlaps(frequbigKaKoar_GRange,lyrGenes,type=c("any"))
frequbigKaKoar_lyrGenes_df=data.frame(as.data.frame(frequbigKaKoar_lyrGenes$frequbigKaKoar_GRange),as.data.frame(frequbigKaKoar_lyrGenes$lyrGenes),test3KaKoar$width)
colnames(frequbigKaKoar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKaKoar_lyrGenes_df<-frequbigKaKoar_lyrGenes_df[(frequbigKaKoar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKaKoar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKaKoar_lyrGenes_df,"KaZar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

Thal_description<-read.xlsx("Lyr_TAIR_Mapman_descript.xlsx")

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

frequbigKaKoar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKaKoar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKaKoar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

UtesMetal=read.xlsx("metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]
Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup=merge(frequbigKaKoar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KaKoar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KaKoar.xlsx",overwrite=T)

Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#119 unique genes
write.xlsx(Thal_MapMan_frequbigKaKoar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaKatoKowaarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Genesarenosacnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

#################################################
#Kowa to Kato arenosa#
#################################################

regionsKKar<-cnvr(resCNMOPSKKar)
singleKKar<-cnvs(resCNMOPSKKar)
testKKar<-findOverlaps(regionsKKar,singleKKar)
test2KKar<- DataFrame(splitAsList(singleKKar$sampleName[subjectHits(testKKar)], queryHits(testKKar)),splitAsList(singleKKar$median[subjectHits(testKKar)], queryHits(testKKar)),splitAsList(singleKKar$mean[subjectHits(testKKar)], queryHits(testKKar)),splitAsList(singleKKar$CN[subjectHits(testKKar)], queryHits(testKKar)))
colnames(test2KKar)<-c("sampleNames","mean","median","CN")
test7KKar<-apply(as.data.frame(test2KKar),1,dupl_cov)
mcols(regionsKKar)<-as.data.frame(as.matrix(test7KKar))

dataKKar<-as.data.frame(regionsKKar)

adddataKKar<-do.call("rbind",apply(dataKKar,1,data_change))
dataKKar<-dataKKar[,-6]
frequKKar<-cbind(dataKKar,adddataKKar)

names(frequKKar)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KKar<-frequKKar[count.fields(textConnection(frequKKar$sampleNames),sep=",")>=5,]

require(openxlsx)
write.xlsx(frequ3KKar,"CnvregionswithmeanarenosaKowaKatoarenosa_minRead6_WL350_minL6_min5samples.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)
frequbigKKar<-read.xlsx("CnvregionswithmeanarenosaKowaKatoarenosa_minRead6_WL350_minL6_min5samples.xlsx",1)
frequbigKKar[,2]<-as.numeric(as.character(frequbigKKar[,2]))
frequbigKKar[,3]<-as.numeric(as.character(frequbigKKar[,3]))

####################################################
#find overlaps between windows and A. lyrata genes#
###################################################
### Convert list of significant markers to a "Genomic Ranges" table. 
frequbigKKar_GRange<-GRanges(seqnames=tolower(frequbigKKar$Scaffold),ranges=IRanges(start=frequbigKKar$start,end=frequbigKKar$end))
values(frequbigKKar_GRange)<-frequbigKKar[,6:11]

testKKar<-findOverlaps(frequbigKKar_GRange,lyrGenes)
test2KKar<-pintersect(frequbigKKar_GRange[queryHits(testKKar)],lyrGenes[subjectHits(testKKar)])
test3KKar<-as.data.frame(test2KKar)
frequbigKKar_lyrGenes=mergeByOverlaps(frequbigKKar_GRange,lyrGenes,type=c("any"))
frequbigKKar_lyrGenes_df=data.frame(as.data.frame(frequbigKKar_lyrGenes$frequbigKKar_GRange),as.data.frame(frequbigKKar_lyrGenes$lyrGenes),test3KKar$width)
colnames(frequbigKKar_lyrGenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKKar_lyrGenes_df<-frequbigKKar_lyrGenes_df[(frequbigKKar_lyrGenes_df$Widthofoverlap>=(0.9*frequbigKKar_lyrGenes_df$gene_size)),]
write.xlsx(frequbigKKar_lyrGenes_df,"ZKaar_temp.xlsx",col.names=TRUE,row.names=FALSE,overwrite=T)

frequbigKKar_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(frequbigKKar_lyrGenes_df,Thal_description, by.x="Gene", by.y="Lyr_Gene",all.x=TRUE)
colnames(frequbigKKar_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup=merge(frequbigKKar_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup<-Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup[,c(2:19,1,20:40)]

write.table(Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KoKaar.txt", sep="\t", row.names=F)
write.xlsx(Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KoKaar.xlsx",overwrite=T)

Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport<-as.numeric(as.character(Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport))
Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup_highsupport<-Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup[Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup$Percentagesupport>=0.9,]
length(unique(Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup_highsupport$Lyr_Gene))
#133 unique genes
write.xlsx(Thal_MapMan_frequbigKKar_lyrGenes_df_desc_w_orthogroup_highsupport,"GenesarenosaKowaKatoarenosacnvregions_minRead6_WL350_minL6_minoverlap90percentofgene_atleast90percentsupport_atleast5samples.xlsx",sheetName="Genesarenosacnmops",col.names=TRUE,row.names=FALSE,overwrite=T)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

KaKoar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KaKoar.xlsx",1)
KKar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KoKaar.xlsx",1)

KaKoar_unique<-KaKoar[!duplicated(KaKoar$Lyr_Gene,KaKoar$Copy_start,KaKoar$Copy_end),]
KKar_unique<-KKar[!duplicated(KKar$Lyr_Gene,KKar$Copy_start,KKar$Copy_end),]

KaKoar_unique$Lyr_Gene<-as.character(KaKoar_unique$Lyr_Gene)
KKar_unique$Lyr_Gene<-as.character(KKar_unique$Lyr_Gene)
KaKoar_unique$CN_class<-as.character(KaKoar_unique$CN_class)
KKar_unique$CN_class<-as.character(KKar_unique$CN_class)

KaKoaronly<-KaKoar_unique[!((KaKoar_unique$Lyr_Gene%in%KKar_unique$Lyr_Gene)&(KaKoar_unique$CN_class==KKar_unique$CN_class[match(KaKoar_unique$Lyr_Gene,KKar_unique$Lyr_Gene)])),]
KKaronly<-KKar_unique[!((KKar_unique$Lyr_Gene%in%KaKoar_unique$Lyr_Gene)&(KKar_unique$CN_class==KaKoar_unique$CN_class[match(KKar_unique$Lyr_Gene,KaKoar_unique$Lyr_Gene)])),]

write.table(KaKoaronly,"KaKoar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(KKaronly,"KoKaar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











