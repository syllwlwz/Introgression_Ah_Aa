options(java.parameters = "-Xmx12000m")
require(openxlsx)

#filelist1 = list.files(path = "./Counts/",pattern = "*.txt",full.names=T) # Make a file list from all count text files
require(gtools)
#filelist1 <- mixedsort(filelist1) # Sort the numbers: 1 to 72

#datafr1 = do.call(cbind,lapply(filelist1,function(fn)read.table(fn,header=FALSE, sep="\t")[,2])) #merge all count files

filelistFC = list.files(path = "./Counts_FC/",pattern = "*.txt$",full.names=T) # Make a file list from all count text files
filelistFC <- mixedsort(filelistFC) # Sort the numbers: 1 to 72

datafrFC = do.call(cbind,lapply(filelistFC,function(fn)read.table(fn,header=TRUE, sep="\t")[,7])) #merge all count files

#sum1<-apply(datafr1,2,"sum")
#sumFC<-apply(datafrFC,2,"sum")
#comp<-data.frame(sum1,sumFC)
#colnames(comp)<-c("Sum_Qualimap_primary_exons","Sum_FeatureCounts_all_exons")
#write.xlsx(comp,"Sum_of_counts_comp_Qualimap_FeatureCounts.xlsx",row.names=F,overwrite=T)

genes <- read.delim(filelistFC[1],header = TRUE,comment.char ="#")[,1] # for the first column
genes2 <- gsub(".v2.1","",as.character(genes))

count0 = cbind(genes2, datafrFC) #merge with gene name
countdata <- count0[,-1]
rownames(countdata) <- count0[,1] 
colnames(countdata) <- gsub(".counts.txt","",gsub("./Counts_FC/","",filelistFC)) # add column names
mode(countdata)
mode(countdata) = "numeric" #convert the class of a matrix from character to numeric
class(countdata)
countdata2<-round(countdata,0)

statsPerSample <- data.frame(t(apply(countdata, 2, summary)))
head(statsPerSample)

setwd("E:/Introgression/Introrna/DESeq2")
# Load the pre-made coldata to represent experimental designs
coldata <- read.delim("Introrna_sample_info.csv", header = TRUE,sep=";")
rownames(coldata) <- coldata$Sample

countdata_df<-data.frame(row.names(countdata),as.data.frame(countdata))
colnames(countdata_df)[1]<-"ID"
#Make a data frame form a matrix
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require(DESeq2)
coldata2<-coldata[match(colnames(countdata),coldata$Sample),]

coldata2$SpeciesPopcomb<-paste0(coldata2$Species,coldata2$Population)
countdata_df2<-countdata_df[coldata2$Treatment=="c",]
countdata3<-countdata2[,coldata2$Treatment=="c"]
coldata3<-coldata2[coldata2$Treatment=="c",]

dds2 <- DESeqDataSetFromMatrix(countData=countdata3,colData=coldata3,design=~Tissue+SpeciesPopcomb)

#Note: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the
#formula and make sure the control level is the first level.

#converting counts to integer mode
#Warnmeldung:
#In DESeqDataSet(se, design = design, ignoreRank) :
#  some variables in design formula are characters, converting to factors

## Check the general properties of the DESeq dataset
print(dds2)
dds_root<-dds2[,dds2$Tissue=="R"]
dds_shoot<-dds2[,dds2$Tissue=="S"]
dds_root$Tissue<-droplevels(dds_root$Tissue)
dds_shoot$Tissue<-droplevels(dds_shoot$Tissue)

#relevel, just specifying the reference level:
dds2$Population <- relevel(as.factor(dds2$Population), ref = "Zapa")

colData(dds_root)<-colData(dds_root)[,c(1:7)]
colData(dds_shoot)<-colData(dds_shoot)[,c(1:7)]

#### Calculate Differential Expression
alpha <- 0.05 #significance level for adjusted p-value, default is 0.1

######################################################################
#Differences between Populations and Species independent of treatment#
######################################################################

design(dds_root)<-~ SpeciesPopcomb
design(dds_shoot)<-~ SpeciesPopcomb
dds_root_calc = DESeq(dds_root,test="Wald",fitType="parametric")
dds_shoot_calc = DESeq(dds_shoot,test="Wald",fitType="parametric")
resultsNames(dds_root_calc)
resultsNames(dds_shoot_calc)

IA_R_thalianaCol0_vs_arenosaMias<-results(dds_root_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaMias"),alpha=0.05)
IA_S_thalianaCol0_vs_arenosaMias<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaMias"),alpha=0.05)
IA_R_thalianaCol0_vs_arenosaZapa<-results(dds_root_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaZapa"),alpha=0.05)
IA_S_thalianaCol0_vs_arenosaZapa<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaZapa"),alpha=0.05)
IA_R_thalianaCol0_vs_halleriMias<-results(dds_root_calc,contrast=c("SpeciesPopcomb","thalianaCol0","halleriMias"),alpha=0.05)
IA_S_thalianaCol0_vs_halleriMias<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","thalianaCol0","halleriMias"),alpha=0.05)
IA_R_thalianaCol0_vs_halleriZapa<-results(dds_root_calc,contrast=c("SpeciesPopcomb","thalianaCol0","halleriZapa"),alpha=0.05)
IA_S_thalianaCol0_vs_halleriZapa<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","thalianaCol0","halleriZapa"),alpha=0.05)

IA_R_thalianaCol0_vs_arenosaMias_sig<-as.data.frame(IA_R_thalianaCol0_vs_arenosaMias[which(IA_R_thalianaCol0_vs_arenosaMias$padj<=0.05&(IA_R_thalianaCol0_vs_arenosaMias$log2FoldChange>=1|IA_R_thalianaCol0_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_S_thalianaCol0_vs_arenosaMias_sig<-as.data.frame(IA_S_thalianaCol0_vs_arenosaMias[which(IA_S_thalianaCol0_vs_arenosaMias$padj<=0.05&(IA_S_thalianaCol0_vs_arenosaMias$log2FoldChange>=1|IA_S_thalianaCol0_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_R_thalianaCol0_vs_arenosaZapa_sig<-as.data.frame(IA_R_thalianaCol0_vs_arenosaZapa[which(IA_R_thalianaCol0_vs_arenosaZapa$padj<=0.05&(IA_R_thalianaCol0_vs_arenosaZapa$log2FoldChange>=1|IA_R_thalianaCol0_vs_arenosaZapa$log2FoldChange<=(-1))),])
IA_S_thalianaCol0_vs_arenosaZapa_sig<-as.data.frame(IA_S_thalianaCol0_vs_arenosaZapa[which(IA_S_thalianaCol0_vs_arenosaZapa$padj<=0.05&(IA_S_thalianaCol0_vs_arenosaZapa$log2FoldChange>=1|IA_S_thalianaCol0_vs_arenosaZapa$log2FoldChange<=(-1))),])
IA_R_thalianaCol0_vs_halleriMias_sig<-as.data.frame(IA_R_thalianaCol0_vs_halleriMias[which(IA_R_thalianaCol0_vs_halleriMias$padj<=0.05&(IA_R_thalianaCol0_vs_halleriMias$log2FoldChange>=1|IA_R_thalianaCol0_vs_halleriMias$log2FoldChange<=(-1))),])
IA_S_thalianaCol0_vs_halleriMias_sig<-as.data.frame(IA_S_thalianaCol0_vs_halleriMias[which(IA_S_thalianaCol0_vs_halleriMias$padj<=0.05&(IA_S_thalianaCol0_vs_halleriMias$log2FoldChange>=1|IA_S_thalianaCol0_vs_halleriMias$log2FoldChange<=(-1))),])
IA_R_thalianaCol0_vs_halleriZapa_sig<-as.data.frame(IA_R_thalianaCol0_vs_halleriZapa[which(IA_R_thalianaCol0_vs_halleriZapa$padj<=0.05&(IA_R_thalianaCol0_vs_halleriZapa$log2FoldChange>=1|IA_R_thalianaCol0_vs_halleriZapa$log2FoldChange<=(-1))),])
IA_S_thalianaCol0_vs_halleriZapa_sig<-as.data.frame(IA_S_thalianaCol0_vs_halleriZapa[which(IA_S_thalianaCol0_vs_halleriZapa$padj<=0.05&(IA_S_thalianaCol0_vs_halleriZapa$log2FoldChange>=1|IA_S_thalianaCol0_vs_halleriZapa$log2FoldChange<=(-1))),])

pdf("Cooks_distance_root_for_outliers.pdf",width=4.72441,height=4.72441,paper="special")
boxplot(log10(assays(dds_root_calc)[["cooks"]]), range=0, las=2, main="Cook's distances") # Check any outliers by Cook's D 
dev.off()
pdf("Cooks_distance_shoot_for_outliers.pdf",width=4.72441,height=4.72441,paper="special")
boxplot(log10(assays(dds_shoot_calc)[["cooks"]]), range=0, las=2, main="Cook's distances") # Check any outliers by Cook's D 
dev.off()

#merge by rownames with thaliana gene description
Thal_description<-read.xlsx("E:/Lyrata_RBH/Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx")
All_lists<-read.xlsx("metal_homeostasis_20200707_v6sorted.xlsx",2)
All_lists<-All_lists[,-1]
All_lists$AGI_Number<-toupper(All_lists$AGI_Number)


###############################
#medofratios#
###############################
#separate root and shoot
#dds_root_norm <- estimateSizeFactors(dds_root)
#normalized_counts_R <-data.frame(counts(dds_root_norm, normalized=TRUE))
#colnames(normalized_counts_R)<-dds_root_norm$Sample.Name

#dds_shoot_norm <- estimateSizeFactors(dds_shoot)
#normalized_counts_S <-data.frame(counts(dds_shoot_norm, normalized=TRUE))
#colnames(normalized_counts_S)<-dds_shoot_norm$Sample.Name

#countdata_df<-cbind(normalized_counts_R,normalized_counts_S)
#Extracted_coldata_R<-as.data.frame(colData(dds_root))
#Extracted_coldata_S<-as.data.frame(colData(dds_shoot))

#colnames(countdata_df)<-c(paste(Extracted_coldata_R$Species,Extracted_coldata_R$Population,Extracted_coldata_R$Treatment,Extracted_coldata_R$Tissue,Extracted_coldata_R$Repeat,sep="_"),paste(Extracted_coldata_S$Species,Extracted_coldata_S$Population,Extracted_coldata_S$Treatment,Extracted_coldata_S$Tissue,Extracted_coldata_S$Repeat,sep="_"))
#countdata_rlog<-countdata_df[,order(colnames(countdata_df))]

#countdata_rlog2<-merge(CDS_df[,c(1,3)],countdata_rlog,by.x="seq_name",by.y="row.names",all.y=T)

#countdata_rlog2[,3:50]<-countdata_rlog2[,3:50]*1000/countdata_rlog2$Tlen
#countdata_rlog_df<-countdata_rlog2[,-2]
#countdata_rlog_df_Tair<-merge(countdata_rlog_df,Thal_description,by.x="seq_name",by.y="Alyr_ID",all.x=TRUE)

#write.table(countdata_rlog_df_Tair,"Countdata_medofratios.table",sep="\t",row.names=F)


#count0 = cbind(genes2, datafrFC) #merge with gene name
#countdata <- count0[,-1]
#rownames(countdata) <- count0[,1] 
#colnames(countdata) <- gsub(".counts.txt","",gsub("./Counts_FC/","",filelistFC)) # add column names
#mode(countdata)
#mode(countdata) = "numeric" #convert the class of a matrix from character to numeric
#class(countdata)

#statsPerSample <- data.frame(t(apply(countdata, 2, summary)))
#head(statsPerSample)

#setwd("E:/Introgression/Introrna/DESeq2")
# Load the pre-made coldata to represent experimental designs
## Calculation of FPKM
# 1) Load gff and Merge
require(Biostrings)
CDS=readDNAStringSet("D:/Lyrata/Alyrata_384_v2.1.cds_primaryTranscriptOnly.fa")
seq_name = names(CDS)
sequence = paste(CDS)
CDS_df <- data.frame(seq_name, sequence)
CDS_df$Tlen<-nchar(as.character(CDS_df$sequence))
test<-sapply(strsplit(as.character(CDS_df$seq_name), "[.]"), `[`, 1)
CDS_df$seq_name<-test
   
FPKM_df = merge(CDS_df[,c(1,3)],countdata2,by.x="seq_name",by.y="row.names",all.y=T)
head(FPKM_df)

# 2) Calculate
  # FPKM=(Counts/(gene size/1000))/(sum mapped fragments/1 million)
for (i in 1:60)
	{Si<-data.frame((FPKM_df[,2+i]/(FPKM_df$Tlen/1000))/(sum(FPKM_df[,2+i])/1000000))
	colnames(Si)<-paste0("FPKM_",colnames(FPKM_df)[2+i])
	FPKM_df<-data.frame(FPKM_df,Si)
	}
colnames(FPKM_df)[1:2]<-c("Alyr_ID","Primary_transcript_length")

#write.table(as.data.frame(FPKM_df),"FPKM_df_named_Ctrl.txt",sep="\t",row.names=FALSE)
#FPKM_df<-read.table("FPKM_df.txt",sep="\t",header=T)

coldata_match<-coldata2[match(coldata2$Sample,colnames(FPKM_df)[3:62]),]
colnames(FPKM_df)[3:62]<-c(paste(coldata_match$Species,coldata_match$Population,coldata_match$Treatment,coldata_match$Tissue,coldata_match$Repeat,sep="_"))
colnames(FPKM_df)[63:122]<-paste("FPKM",colnames(FPKM_df)[3:62],sep="_")

countdata_FPKM<-data.frame(FPKM_df[,1:2],FPKM_df[,3:62][,order(colnames(FPKM_df[3:62]))],FPKM_df[,63:122][,order(colnames(FPKM_df[63:122]))])

#add gff description and info
gff<-read.table("D:/Lyrata/Alyrata_384_v2.1.gene.gff3",sep="\t",header=F)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
gff<-gff[gff$Type=="gene",]
ID<-substr(as.character(gff$Attributes),4,nchar(as.character(gff$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
gff$Lyr_Gene<-ID3[,1]

IA_R_thalianaCol0_vs_arenosaMias_sig_lyr<-merge(IA_R_thalianaCol0_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_thalianaCol0_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_R_thalianaCol0_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted<-IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_R_thalianaCol0_vs_arenosaMias_inclAth_Ctrl.table",sep="\t",row.names=F)


IA_S_thalianaCol0_vs_arenosaMias_sig_lyr<-merge(IA_S_thalianaCol0_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_thalianaCol0_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_S_thalianaCol0_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted<-IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_S_thalianaCol0_vs_arenosaMias_inclAth_Ctrl.table",sep="\t",row.names=F)


IA_R_thalianaCol0_vs_arenosaZapa_sig_lyr<-merge(IA_R_thalianaCol0_vs_arenosaZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_thalianaCol0_vs_arenosaZapa_sig_lyr)[1]<-"Alyr_ID"
IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc<-merge(Thal_description,IA_R_thalianaCol0_vs_arenosaZapa_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted<-IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc[order(-abs(IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc$log2FoldChange)),]

IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff<-merge(IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff<-IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff_lists,"IA_R_thalianaCol0_vs_arenosaZapa_inclAth_Ctrl.table",sep="\t",row.names=F)


IA_S_thalianaCol0_vs_arenosaZapa_sig_lyr<-merge(IA_S_thalianaCol0_vs_arenosaZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_thalianaCol0_vs_arenosaZapa_sig_lyr)[1]<-"Alyr_ID"
IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc<-merge(Thal_description,IA_S_thalianaCol0_vs_arenosaZapa_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted<-IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc[order(-abs(IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc$log2FoldChange)),]

IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff<-merge(IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff<-IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff_lists,"IA_S_thalianaCol0_vs_arenosaZapa_inclAth_Ctrl.table",sep="\t",row.names=F)


IA_R_thalianaCol0_vs_halleriMias_sig_lyr<-merge(IA_R_thalianaCol0_vs_halleriMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_thalianaCol0_vs_halleriMias_sig_lyr)[1]<-"Alyr_ID"
IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc<-merge(Thal_description,IA_R_thalianaCol0_vs_halleriMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted<-IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc[order(-abs(IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc$log2FoldChange)),]

IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff<-merge(IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff<-IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff_lists,"IA_R_thalianaCol0_vs_halleriMias_inclAth_Ctrl.table",sep="\t",row.names=F)


IA_S_thalianaCol0_vs_halleriMias_sig_lyr<-merge(IA_S_thalianaCol0_vs_halleriMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_thalianaCol0_vs_halleriMias_sig_lyr)[1]<-"Alyr_ID"
IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc<-merge(Thal_description,IA_S_thalianaCol0_vs_halleriMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted<-IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc[order(-abs(IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc$log2FoldChange)),]

IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff<-merge(IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff<-IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_thalianaCol0_vs_halleriMias_sig_Ath_desc_sorted_gff_lists,"IA_S_thalianaCol0_vs_halleriMias_inclAth_Ctrl.table",sep="\t",row.names=F)


IA_R_thalianaCol0_vs_halleriZapa_sig_lyr<-merge(IA_R_thalianaCol0_vs_halleriZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_thalianaCol0_vs_halleriZapa_sig_lyr)[1]<-"Alyr_ID"
IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc<-merge(Thal_description,IA_R_thalianaCol0_vs_halleriZapa_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted<-IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc[order(-abs(IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc$log2FoldChange)),]

IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff<-merge(IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff<-IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff_lists,"IA_R_thalianaCol0_vs_halleriZapa_inclAth_Ctrl.table",sep="\t",row.names=F)


IA_S_thalianaCol0_vs_halleriZapa_sig_lyr<-merge(IA_S_thalianaCol0_vs_halleriZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_thalianaCol0_vs_halleriZapa_sig_lyr)[1]<-"Alyr_ID"
IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc<-merge(Thal_description,IA_S_thalianaCol0_vs_halleriZapa_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted<-IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc[order(-abs(IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc$log2FoldChange)),]

IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff<-merge(IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff<-IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_thalianaCol0_vs_halleriZapa_sig_Ath_desc_sorted_gff_lists,"IA_S_thalianaCol0_vs_halleriZapa_inclAth_Ctrl.table",sep="\t",row.names=F)








#compare

new_R<-read.xlsx("Deseq2_introgression_inclAth_Ctrl.xlsx",1)
new_S<-read.xlsx("Deseq2_introgression_inclAth_Ctrl.xlsx",2)

OldZ_R<-read.xlsx("Deseq2_introgression_MAa_vs_ZAa_wo_MAa_vs_ZAh_inclAth.xlsx",1)
OldZ_S<-read.xlsx("Deseq2_introgression_MAa_vs_ZAa_wo_MAa_vs_ZAh_inclAth.xlsx",2)

OldM_R<-read.xlsx("Deseq2_introgression_inclAth.xlsx",1)
OldM_S<-read.xlsx("Deseq2_introgression_inclAth.xlsx",2)

Old_or_R_comb<-rbind(OldZ_R,OldM_R)
Old_or_R<-Old_or_R_comb[!duplicated(Old_or_R_comb),]
#2523
Old_or_S_comb<-rbind(OldZ_S,OldM_S)
Old_or_S<-Old_or_S_comb[!duplicated(Old_or_S_comb),]
#1709

length(new_R$Alyr_ID[!duplicated(new_R$Alyr_ID)])
#2524
length(new_R$Name[!duplicated(new_R$Name)])
#50
paste0(new_R$Name[!duplicated(new_R$Name)],sep=",",collapse="")
#"NA,FRO2,ZIP5,ZIP4,MTP10,SPL4,CAX10,NAS4,IRT3,PDF1.1,VTL2,MTP9,PDF2.1,MCA2,HMA4,GLR2.3,GLR2.2,COPT6,ZIP6,CHX21,COPT4,MTP11,PDF-like?,MTP1,FRD3,LAC7,MT2a,FER2,ZIP1,F6'H1,CNGC19,CNGC20,VTL5,CAX3,CHX20,OSCA2.5,OSCA1.6,IRT2,HMA2,HMA3,ZIP9,FER1,ZIFl1,FRO4,FRO5,MCW1,BGLU42,MCU5,ZIP8,TIL,"
50/2524*100
#[1] 1.980983

length(OldZ_R$Alyr_ID[!duplicated(OldZ_R$Alyr_ID)])
#1816
length(OldZ_R$Name[!duplicated(OldZ_R$Name)])
#46
paste0(OldZ_R$Name[!duplicated(OldZ_R$Name)],sep=",",collapse="")
#"NA,FRO2,ZIP4,PCR2,CGNC7,CGNC8,LPR1,ZIP10,SPL4,CAX10,OSCA2.1,IRT3,PDF2.1,MCA2,HMA4,GLR2.2,COPT6,ZIP6,CIPK11,MTM1-like,CGNC3,MTP1,OAS-TL,FRD3,LAC7,MT2a,COX17-1,CNGC20,OSCA2.5,bHLH039,OSCA1.6,IRT2,HMA2,FER1,HTR15,ZIFl1,FRO4,MCW1,BGLU42,MCU5,ZIP8,FRO6,FRO7,SPL-like,TIL,ZIP12,"
46/1816*100
#[1] 2.53304

length(OldM_R$Alyr_ID[!duplicated(OldM_R$Alyr_ID)])
#1976
length(OldM_R$Name[!duplicated(OldM_R$Name)])
#42
paste0(OldM_R$Name[!duplicated(OldM_R$Name)],sep=",",collapse="")
#"NA,FRO2,ZIP4,PCR2,VTL1,LPR1,SPL4,CAX10,PDF1.5,IRT3,CHX16,VTL2,MTP9,PDF2.1,MCA2,GLR2.3,COPT6,ZIP6,COPT4,MTP11,MTM1-like,CGNC3,MTP1,OAS-TL,FRD3,LAC7,ZIP1,COX17-1,CNGC20,VTL5,CAX3,OSCA2.5,bHLH039,CAX4,FER1,HTR15,FRO4,MCW1,BGLU42,SPL-like,TIL,ZIP12,"
42/1976*100
#[1] 2.125506

length(Old_or_R_comb$Alyr_ID[!duplicated(Old_or_R_comb$Alyr_ID)])
#2523
length(Old_or_R_comb$Name[!duplicated(Old_or_R_comb$Name)])
#58
paste0(Old_or_R_comb$Name[!duplicated(Old_or_R_comb$Name)],sep=",",collapse="")
#"NA,FRO2,ZIP4,PCR2,CGNC7,CGNC8,LPR1,ZIP10,SPL4,CAX10,OSCA2.1,IRT3,PDF2.1,MCA2,HMA4,GLR2.2,COPT6,ZIP6,CIPK11,MTM1-like,CGNC3,MTP1,OAS-TL,FRD3,LAC7,MT2a,COX17-1,CNGC20,OSCA2.5,bHLH039,OSCA1.6,IRT2,HMA2,FER1,HTR15,ZIFl1,FRO4,MCW1,BGLU42,MCU5,ZIP8,FRO6,FRO7,SPL-like,TIL,ZIP12,VTL1,PDF1.5,CHX16,VTL2,MTP9,GLR2.3,COPT4,MTP11,ZIP1,VTL5,CAX3,CAX4,"
58/2523*100
#[1] 2.298851

length(Cand_R_samedir$Alyr_ID[!duplicated(Cand_R_samedir$Alyr_ID)])
#2065
length(Cand_R_samedir$Name.x[!duplicated(Cand_R_samedir$Name.x)])
#44
paste0(Cand_R_samedir$Name.x[!duplicated(Cand_R_samedir$Name.x)],sep=",",collapse="")
#"NA,FRO2,ZIP4,HMA2,MTP10,SPL4,IRT3,PDF1.1,VTL2,MTP9,FRD3,LAC7,MT2a,FER2,ZIP1,F6'H1,CNGC20,VTL5,MCA2,HMA4,GLR2.3,GLR2.2,COPT6,ZIP6,MTP11,PDF-like?,MTP1,PDF2.1,COPT4,CAX3,OSCA2.5,FER1,ZIFl1,FRO4,FRO5,CAX10,ZIP9,IRT2,OSCA1.6,BGLU42,ZIP8,MCU5,MCW1,TIL,"
44/2065*100
#[1] 2.130751










length(new_S$Alyr_ID[!duplicated(new_S$Alyr_ID)])
#1608
length(new_S$Name[!duplicated(new_S$Name)])
#24
paste0(new_S$Name[!duplicated(new_S$Name)],sep=",",collapse="")
# "NA,FRO2,ZIP5,OSCA1.3,CAX10,IRT3,PDF2.4,CHX2,PDF2.1,HMA4,unknown,ZIP6,CHX21,COPT4,MTP1,FER2,COPT2,OSCA2.5,MTP2,OSCA1.2,HMA2,HMA3,PDF1.2C,ZIP8,"
24/1608*100
#[1] 1.492537

length(OldZ_S$Alyr_ID[!duplicated(OldZ_S$Alyr_ID)])
#1283
length(OldZ_S$Name[!duplicated(OldZ_S$Name)])
#26
paste0(OldZ_S$Name[!duplicated(OldZ_S$Name)],sep=",",collapse="")
#"NA,CNGC10,FRO2,GLR3.4,CCH-like,CAX10,OSCA2.1,IRT3,BTSl1,PDF2.1,HMA4,unknown,ZIP6,MTM1-like,FER2,COPT2,bHLH039,MTP2,OSCA1.2,MTM1,HMA2,HMA3,ZIP9,CAX7,PDF1.2C,ZIP8,"
26/1283*100
#[1] 2.0265

length(OldM_S$Alyr_ID[!duplicated(OldM_S$Alyr_ID)])
#1306
length(OldM_S$Name[!duplicated(OldM_S$Name)])
#24
paste0(OldM_S$Name[!duplicated(OldM_S$Name)],sep=",",collapse="")
#"NA,CNGC10,FRO2,GLR3.4,ZIP4,CAX10,IRT3,PDF2.1,unknown,ZIP6,COPT4,MTM1-like,ZIP1,SPL5,COPT2,OSCA2.5,MTP2,OSCA1.2,MTM1,ZIP9,CAX7,MSL9,PDF1.2C,ZIP8,"
24/1306*100
#[1] 1.837672

length(Old_or_S_comb$Alyr_ID[!duplicated(Old_or_S_comb$Alyr_ID)])
#1709
length(Old_or_S_comb$Name[!duplicated(Old_or_S_comb$Name)])
#32
paste0(Old_or_S_comb$Name[!duplicated(Old_or_S_comb$Name)],sep=",",collapse="")
#"NA,CNGC10,FRO2,GLR3.4,CCH-like,CAX10,OSCA2.1,IRT3,BTSl1,PDF2.1,HMA4,unknown,ZIP6,MTM1-like,FER2,COPT2,bHLH039,MTP2,OSCA1.2,MTM1,HMA2,HMA3,ZIP9,CAX7,PDF1.2C,ZIP8,ZIP4,COPT4,ZIP1,SPL5,OSCA2.5,MSL9,"
32/1709*100
#[1] 1.87244

length(Cand_S_samedir$Alyr_ID[!duplicated(Cand_S_samedir$Alyr_ID)])
#1320
length(Cand_S_samedir$Name.x[!duplicated(Cand_S_samedir$Name.x)])
#19
paste0(Cand_S_samedir$Name.x[!duplicated(Cand_S_samedir$Name.x)],sep=",",collapse="")
#"NA,HMA2,PDF1.2C,PDF2.4,IRT3,HMA4,unknown,ZIP6,MTP1,PDF2.1,COPT4,COPT2,OSCA2.5,MTP2,CAX10,FRO2,HMA3,OSCA1.2,ZIP8,"
19/1320*100
#[1] 1.439394



























Log<-rlog(dds2, blind = FALSE)
LogS <- rlog(dds_shoot, blind = FALSE)
LogR <- rlog(dds_root, blind = FALSE)
pdf("PCA_all_introgression.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA(Log, intgroup=c("Population", "Tissue","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=as.factor(pcaData$Species),pch=c(16,17)[as.factor(pcaData$Tissue)],xlab="",ylab="",main="Total",cex=1.5,cex.lab=1.2,axes=F)
legend("topright",inset=c(-0.22,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Root","Shoot"),pch=c(16,17),title="Tissue",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))
pcaData <- plotPCA(LogR, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Root",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA(LogS, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.5,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.35,0.3),legend=c("Col0","Mias","Zapa"),pch=c(18,16,17),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)

dev.off()


#######################################################
#Other PC axes#
#######################################################
rv <- rowVars(assay(Log))
select <- order(rv, decreasing = TRUE)[seq_len(min(31073,length(rv)))]
pca <- prcomp(t(assay(Log)[select, ]))
ev <- pca$sdev^2
ev[ev> mean(ev)]
length(ev[ev> mean(ev)])
#7

rv <- rowVars(assay(LogR))
select <- order(rv, decreasing = TRUE)[seq_len(min(31073,length(rv)))]
pca <- prcomp(t(assay(LogR)[select, ]))
ev <- pca$sdev^2
ev[ev> mean(ev)]
length(ev[ev> mean(ev)])
#4

rv <- rowVars(assay(LogS))
select <- order(rv, decreasing = TRUE)[seq_len(min(31073,length(rv)))]
pca <- prcomp(t(assay(LogS)[select, ]))
ev <- pca$sdev^2
ev[ev> mean(ev)]
length(ev[ev> mean(ev)])
#4

library(genefilter)
library(ggplot2)
library(ggrepel)

plotPCA.PC2_3 <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3],group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(2:3)]
    return(d)
  }
    ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 

}

plotPCA.PC1_3 <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3],group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(1,3)]
    return(d)
  }
    ggplot(data = d, aes_string(x = "PC1", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 

}

plotPCA.PC3_4 <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4],group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(3,4)]
    return(d)
  }
    ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 

}


pdf("PCA_all_introgression_PC1_3.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC1_3(Log, intgroup=c("Population", "Tissue","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC3~pcaData$PC1,col=as.factor(pcaData$Species),pch=c(16,17)[as.factor(pcaData$Tissue)],xlab="",ylab="",main="Total",cex=1.5,cex.lab=1.2,axes=F)
legend("topright",inset=c(-0.22,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Root","Shoot"),pch=c(16,17),title="Tissue",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))
pcaData <- plotPCA.PC1_3(LogR, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC3~pcaData$PC1,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Root",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC1_3(LogS, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC3~pcaData$PC1,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.5,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.35,0.3),legend=c("Col0","Mias","Zapa"),pch=c(18,16,17),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)

dev.off()


pdf("PCA_all_introgression_PC2_3.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC2_3(Log, intgroup=c("Population", "Tissue","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC3~pcaData$PC2,col=as.factor(pcaData$Species),pch=c(16,17)[as.factor(pcaData$Tissue)],xlab="",ylab="",main="Total",cex=1.5,cex.lab=1.2,axes=F)
legend("topright",inset=c(-0.22,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Root","Shoot"),pch=c(16,17),title="Tissue",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC2 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))
pcaData <- plotPCA.PC2_3(LogR, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC3~pcaData$PC2,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Root",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC2 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC2_3(LogS, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC3~pcaData$PC2,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC2 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.5,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.35,0.3),legend=c("Col0","Mias","Zapa"),pch=c(18,16,17),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)

dev.off()

pdf("PCA_all_introgression_PC3_4.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC3_4(Log, intgroup=c("Population", "Tissue","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC4~pcaData$PC3,col=as.factor(pcaData$Species),pch=c(16,17)[as.factor(pcaData$Tissue)],xlab="",ylab="",main="Total",cex=1.5,cex.lab=1.2,axes=F)
legend("topright",inset=c(-0.22,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Root","Shoot"),pch=c(16,17),title="Tissue",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC3 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC4 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))
pcaData <- plotPCA.PC3_4(LogR, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC4~pcaData$PC3,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Root",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC3 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC4 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC3_4(LogS, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC4~pcaData$PC3,col=as.factor(pcaData$Species),pch=c(18,16,17)[as.factor(pcaData$Population)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC3 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC4 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.5,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,fill=c("red","black","green"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.35,0.3),legend=c("Col0","Mias","Zapa"),pch=c(18,16,17),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)

dev.off()

################################
#Subset#
################################
LogSMZ <-LogS[,LogS$Species=="arenosa"|LogS$Species=="halleri"]
LogRMZ <-LogR[,LogR$Species=="arenosa"|LogR$Species=="halleri"]


ev <- pca$sdev^2
ev[ev> mean(ev)]
#4

library(genefilter)
library(ggplot2)
library(ggrepel)

plotPCA.PC2_3 <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3],group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(2:3)]
    return(d)
  }
    ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 

}

plotPCA.PC1_3 <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3],group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(1,3)]
    return(d)
  }
    ggplot(data = d, aes_string(x = "PC1", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 

}

pdf("PCA_R_MZ_introgression.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA(LogSMZ, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,bg=c("green","blue")[as.factor(pcaData$Species)],pch=c(21,22)[as.factor(pcaData$Population)],col=c("black","red")[as.factor(pcaData$Treatment)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F,lwd=2)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.22,0),legend=c("A. arenosa","A. halleri"),text.font=3,fill=c("green","blue"),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Mias","Zapa"),pch=c(21,22),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
legend("topright",inset=c(-0.22,0.8),legend=c("Control","Treatment"),col=c("black","red"),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0,pch=22)
par(mar=c(3,3,2,1))

pcaData13 <- plotPCA.PC1_3(LogSMZ, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
plot(pcaData13$PC3~pcaData13$PC1,bg=c("green","blue")[as.factor(pcaData13$Species)],pch=c(21,22)[as.factor(pcaData13$Population)],col=c("black","red")[as.factor(pcaData13$Treatment)],xlab="",ylab="",main="",cex=1.5,cex.lab=1.2,axes=F,lwd=2)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
percentVar <- round(100 * attr(pcaData13, "percentVar"))
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))

pcaData23 <- plotPCA.PC2_3(LogSMZ, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
plot(pcaData23$PC3~pcaData23$PC2,bg=c("green","blue")[as.factor(pcaData23$Species)],pch=c(21,22)[as.factor(pcaData23$Population)],col=c("black","red")[as.factor(pcaData23$Treatment)],xlab="",ylab="",main="",cex=1.5,cex.lab=1.2,axes=F,lwd=2)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
percentVar <- round(100 * attr(pcaData23, "percentVar"))
mtext(side=1,line=2,paste0("PC2 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC3 (",percentVar[2],"%) "))
box()

dev.off()


###################################################
#Hierarchical clustering#
###################################################

#Calculate sample distance
sampleDists <- dist( t( assay(LogS) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(dds_shoot$Species,dds_shoot$Population,dds_shoot$Treatment,dds_shoot$Repeat,sep="-" )
colnames(sampleDistMatrix) <- NULL

attributes(sampleDists)$Labels<-paste(dds_shoot$Species,dds_shoot$Population,dds_shoot$Treatment,dds_shoot$Repeat, sep="-" )
hcluster<-hclust(sampleDists,method="ward.D") 
pdf("Hierarchical_clustering_introgression_rlog_shoot.pdf",width=8,height=8,paper="special")
plot(hcluster) 
dev.off()

require("RColorBrewer")
require("gplots")

pdf("Hierarchical_clustering_with_heatmap_introgression_rlog_shoot.pdf",width=8,height=8,paper="special")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins = c(5,10))
dev.off()

sampleDists <- dist( t( assay(LogR) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(dds_root$Species,dds_root$Population,dds_root$Treatment,dds_root$Repeat, sep="-" )
colnames(sampleDistMatrix) <- NULL

attributes(sampleDists)$Labels<-paste(dds_root$Species,dds_root$Population,dds_root$Treatment,dds_root$Repeat, sep="-" )
hcluster<-hclust(sampleDists,method="ward.D") 
pdf("Hierarchical_clustering_introgression_rlog_root.pdf",width=8,height=8,paper="special")
plot(hcluster) 
dev.off()
pdf("Hierarchical_clustering_with_heatmap_introgression_rlog_root.pdf",width=8,height=8,paper="special")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins = c(5,10))
dev.off()


sampleDists <- dist( t( assay(Log) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(dds2$Species,dds2$Population,dds2$Treatment,dds2$Repeat,dds2$Tissue,sep="-" )
colnames(sampleDistMatrix) <- NULL

attributes(sampleDists)$Labels<-paste(dds2$Species,dds2$Population,dds2$Treatment,dds2$Repeat,dds2$Tissue, sep="-" )
hcluster<-hclust(sampleDists,method="ward.D") 
pdf("Hierarchical_clustering_introgression_rlog_total.pdf",width=4.72441,height=4.72441,paper="special")
plot(hcluster) 
dev.off()
pdf("Hierarchical_clustering_with_heatmap_introgression_rlog_total.pdf",width=4.72441,height=4.72441,paper="special")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins = c(5,10))
dev.off()

library("RColorBrewer")
library("pheatmap")
LogR <- rlog(dds_root, blind = FALSE)
LogS <- rlog(dds_shoot, blind = FALSE)







sampleDistsR <- dist(t(assay(LogR)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- Log$Sample.Name
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Hierarchical_clustering_Cu_rlog.pdf",width=4.72441,height=4.72441,paper="special")
pheatmap(as.matrix(data.frame()),
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
#biocLite("genefilter")
library( "genefilter" )
library( "gplots" )
library( "RColorBrewer" )

# Heatmaps are more comparable if the expression data has been transformed. We can transform our data using the function rlog()
rld <- rlog( dds_root )
colnames(rld)<-dds_root$Sample.Name
# We then generate a heatmap by assaying the rld data based on our adjusted p values. Other factors could also be used to subet the gene list.

Krom_Cu_Ah<-merge(IA_R_Krom_sig,countdata_rlog_df,by.x="row.names",by.y="Alyr_Gene",all.x=TRUE)
Krom_R_Cu_Ah<-Krom_Cu_Ah[,c(1:19)]
Krom_R_Cu_Ah_ordered<-Krom_R_Cu_Ah[,c("Row.names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Krom_Ctrl_R1","Krom_Ctrl_R2","Krom_Ctrl_R3","Krom_Cu_R1","Krom_Cu_R2","Krom_Cu_R3","Kosi_Ctrl_R1","Kosi_Ctrl_R2","Kosi_Ctrl_R3","Kosi_Cu_R1","Kosi_Cu_R2","Kosi_Cu_R3")]
colnames(Krom_R_Cu_Ah_ordered)[8:19]<-c("Krom control R1","Krom control R2","Krom control R3","Krom high Cu R1","Krom high Cu R2","Krom high Cu R3","Kosi control R1","Kosi control R2","Kosi control R3","Kosi high Cu R1","Kosi high Cu R2","Kosi high Cu R3")
pheatmap(Krom_R_Cu_Ah_ordered[,8:19], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(bold(Krom~root~Cu~reaction)))


Intro1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists,Kosi_R_Cu_Ath_desc_sorted_gff_lists)
Intro<-Intro1[!duplicated(Intro1$Lyr_Gene),]
pheatmap(Intro[,26:37], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(bold(Cu~reaction)))

IntroS<-Kosi_S_Cu_Ath_desc_sorted_gff_lists[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists$Lyr_Gene),]
pheatmap(IntroS[,38:49], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(bold(Cu~reaction)),treeheight_row = 0, treeheight_col = 0)




Intro1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists,Kosi_R_Cu_Ath_desc_sorted_gff_lists)
Intro<-Intro1[!duplicated(Intro1$Lyr_Gene),]
rownames(Intro)<-c("COPT2","FSD1","FRO4_1","FRO4_2","CITF1","CCH","UNE15","YSL2","ZIP2")
R<-pheatmap(Intro[,c(35:37,32:34,29:31,26:28)], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Root),treeheight_ro =0,treeheight_col=0,gaps_col=c(3,6,9),show_colnames=F,fontsize_row=10,fontsize_col=10,cellwidth=35,cellheight=35)$gtable

IntroS<-Kosi_S_Cu_Ath_desc_sorted_gff_lists[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists$Lyr_Gene),]
rownames(IntroS)<-c("CCS","CV","CSD2","PKP4")
labels_col<-c("","Kosi -Cu","","","Kosi +Cu","","","Krom -Cu","","","Krom +Cu","")
S<-pheatmap(IntroS[,c(47:49,44:46,41:43,38:40)], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Shoot),treeheight_ro =0,treeheight_col=0,gaps_col=c(3,6,9),labels_col=labels_col,angle_col=0,fontsize_row=10,fontsize_col=10,cellwidth=35,cellheight=35)$gtable

require(gridExtra)
library(gtable)
library(grid)
g <- rbind(R, S,size="first")
g$widths <- unit.pmax(R$widths, S$widths)
grid.newpage()
pdf("Heatmap_Cu_reaction_stricter.pdf",width=9.44882,height=9.44882,paper="special")
grid.draw(g)
dev.off()


#mean and side by side
Intro1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists,Kosi_R_Cu_Ath_desc_sorted_gff_lists)
Intro2<-Intro1[!duplicated(Intro1$Lyr_Gene),]
Intro_KosiCtrl<-apply(Intro2[,c(35:37)],1,mean)
Intro_KosiCu<-apply(Intro2[,c(32:34)],1,mean)
Intro_KromCtrl<-apply(Intro2[,c(29:31)],1,mean)
Intro_KromCu<-apply(Intro2[,c(26:28)],1,mean)
Intro<-cbind(Intro_KosiCtrl,Intro_KosiCu,Intro_KromCtrl,Intro_KromCu)
rownames(Intro)<-c("COPT2","FSD1","FRO4","FRO5","CITF1","CCH","UNE15","YSL2","ZIP2")
labels_col<-c("",expression(paste("Kosi ",bold(L),"Cu")),"","",expression(paste("Kosi ",bold(H),"Cu")),"","",expression(paste("Krom ",bold(L),"Cu")),"","",expression(paste("Krom ",bold(H),"Cu")),"")
labels_row<-c(expression(italic(COPT2^"*")),expression(italic(FSD1^"*")),expression(italic(FRO4^"*")),expression(italic(FRO5^"*")),expression(italic(CITF1)),expression(italic(CCH)),expression(italic(UNE15)),expression(italic(YSL2)),expression(italic(ZIP2)))

R<-pheatmap(Intro, scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Root),treeheight_row =0,treeheight_col=0,show_colnames=T,angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,labels_col=labels_col,labels_row=labels_row,fontsize=12)$gtable

IntroS1<-Kosi_S_Cu_Ath_desc_sorted_gff_lists[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists$Lyr_Gene),]
IntroS_KosiCtrl<-apply(IntroS1[,c(47:49)],1,mean)
IntroS_KosiCu<-apply(IntroS1[,c(44:46)],1,mean)
IntroS_KromCtrl<-apply(IntroS1[,c(41:43)],1,mean)
IntroS_KromCu<-apply(IntroS1[,c(38:40)],1,mean)
IntroS<-cbind(IntroS_KosiCtrl,IntroS_KosiCu,IntroS_KromCtrl,IntroS_KromCu)
rownames(IntroS)<-c("CCS","CV","CSD2","PKP4")
labels_row<-c(expression(italic(CCS)),expression(italic(CV)),expression(italic(CSD2)),expression(italic(PKP4)))
labels_col<-c("",expression(paste("Kosi ",bold(L),"Cu")),"","",expression(paste("Kosi ",bold(H),"Cu")),"","",expression(paste("Krom ",bold(L),"Cu")),"","",expression(paste("Krom ",bold(H),"Cu")),"")
IntroS<-IntroS[c(1,3,2,4),]
S<-pheatmap(IntroS, scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(Shoot),treeheight_row=0,treeheight_col=0,labels_col=labels_col,angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,labels_row=labels_row,fontsize=12)$gtable

require(gridExtra)
library(gtable)
library(grid)
g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Cu_reaction_stricter_beside_smaller.pdf",width=6,height=6,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()

























Intro1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists,Kosi_R_Cu_Ath_desc_sorted_gff_lists)
Intro<-Intro1[!duplicated(Intro1$Lyr_Gene),]
rownames(Intro)<-c("COPT2","FSD1","FRO4_1","FRO4_2","CITF1","CCH","UNE15","YSL2","ZIP2")
R<-pheatmap(Intro[,c(35:37,32:34,29:31,26:28)], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Root),treeheight_ro =0,treeheight_col=0,gaps_col=c(3,6,9),show_colnames=F,fontsize_row=10,fontsize_col=10,cellwidth=35,cellheight=35)$gtable

IntroS<-Kosi_S_Cu_Ath_desc_sorted_gff_lists[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists$Lyr_Gene),]
rownames(IntroS)<-c("CCS","CV","CSD2","PKP4")
labels_col<-c("","Kosi -Cu","","","Kosi +Cu","","","Krom -Cu","","","Krom +Cu","")
S<-pheatmap(IntroS[,c(47:49,44:46,41:43,38:40)], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Shoot),treeheight_ro =0,treeheight_col=0,gaps_col=c(3,6,9),labels_col=labels_col,angle_col=0,fontsize_row=10,fontsize_col=10,cellwidth=35,cellheight=35)$gtable

require(gridExtra)
library(gtable)
library(grid)
g <- rbind(R, S,size="first")
g$widths <- unit.pmax(R$widths, S$widths)
grid.newpage()
pdf("Heatmap_Cu_reaction_stricter.pdf",width=9.44882,height=9.44882,paper="special")
grid.draw(g)
dev.off()


#mean and side by side
Intro1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists,Kosi_R_Cu_Ath_desc_sorted_gff_lists)
Intro2<-Intro1[!duplicated(Intro1$Lyr_Gene),]
Intro_KosiCtrl<-apply(Intro2[,c(35:37)],1,mean)
Intro_KosiCu<-apply(Intro2[,c(32:34)],1,mean)
Intro_KromCtrl<-apply(Intro2[,c(29:31)],1,mean)
Intro_KromCu<-apply(Intro2[,c(26:28)],1,mean)
Intro<-cbind(Intro_KosiCtrl,Intro_KosiCu,Intro_KromCtrl,Intro_KromCu)
rownames(Intro)<-c("COPT2","FSD1","FRO4","FRO5","CITF1","CCH","UNE15","YSL2","ZIP2")
labels_col<-c("",expression(paste("Kosi ",bold(L),"Cu")),"","",expression(paste("Kosi ",bold(H),"Cu")),"","",expression(paste("Krom ",bold(L),"Cu")),"","",expression(paste("Krom ",bold(H),"Cu")),"")
labels_row<-c(expression(italic(COPT2^"*")),expression(italic(FSD1^"*")),expression(italic(FRO4^"*")),expression(italic(FRO5^"*")),expression(italic(CITF1)),expression(italic(CCH)),expression(italic(UNE15)),expression(italic(YSL2)),expression(italic(ZIP2)))





library("RColorBrewer")
library("pheatmap")
library( "genefilter" )
library( "gplots" )

FPKM_df<-read.table("FPKM_df.txt",sep="\t",header=T)

colnames_heatmap1<-coldata2[match(gsub("FPKM_","",colnames(FPKM_df)[63:122]),coldata2$Sample),c(3,5,6,2)]
colnames_heatmap<-apply(colnames_heatmap1,1,paste, collapse="-")
#colnames(FPKM_df)[3:62]<-colnames_heatmap
#colnames(FPKM_df)[63:122]<-colnames_heatmap
FPKM2<-FPKM_df[,63:122]
colnames(FPKM2)<-colnames_heatmap

FPKM_mean<-t(apply(FPKM2,1,function(x) tapply(x,colnames(FPKM2),mean)))
rownames(FPKM_mean)<-FPKM_df$ Alyr_ID

#AL3G52820 HMA4
#AL4G46270 MTP1
#AL1G54170 HMA2
#AL3G19430 FRD3
#AL6G35300 FRO4
#AL4G24850 ZIP6
#AL7G22080 HMA3

#"HMA2","FRD3","HMA4","ZIP6","MTP1","FRO4","HMA3"
require(gridExtra)
library(gtable)
library(grid)

R<-pheatmap(FPKM_mean[rownames(FPKM_mean)=="AL3G52820"|rownames(FPKM_mean)=="AL4G46270"|rownames(FPKM_mean)=="AL1G54170"|rownames(FPKM_mean)=="AL3G19430"|rownames(FPKM_mean)=="AL4G24850"|rownames(FPKM_mean)=="AL7G22080",1:10],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Root",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col=colnames(FPKM_mean)[1:10],angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,fontsize=12,
labels_row=c(expression(italic("HMA2")),expression(italic("FRD3")),expression(italic("HMA4")),expression(italic("ZIP6")),expression(italic("MTP1")),expression(italic("HMA3"))))$gtable
S<-pheatmap(FPKM_mean[rownames(FPKM_mean)=="AL3G52820"|rownames(FPKM_mean)=="AL4G46270"|rownames(FPKM_mean)=="AL1G54170"|rownames(FPKM_mean)=="AL3G19430"|rownames(FPKM_mean)=="AL4G24850"|rownames(FPKM_mean)=="AL7G22080",11:20],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Shoot",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col=colnames(FPKM_mean)[11:20],angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,fontsize=12,
labels_row=c(expression(italic("HMA2")),expression(italic("FRD3")),expression(italic("HMA4")),expression(italic("ZIP6")),expression(italic("MTP1")),expression(italic("HMA3"))))$gtable

g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Introgression_HM_cand.pdf",width=12,height=12,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()



#AL1G10810 FRO2
#AL5G25930 COPT2
#AL5G38430 bHLH39
#AL6G35300 FRO4
#AL7G29170 YSL1
#AL7G34350 IRT1

#"FRO2","COPT2","bHLH39","FRO4","YSL1","IRT1"

R<-pheatmap(FPKM_mean[rownames(FPKM_mean)=="AL1G10810"|rownames(FPKM_mean)=="AL5G25930"|rownames(FPKM_mean)=="AL5G38430"|rownames(FPKM_mean)=="AL6G35300"|rownames(FPKM_mean)=="AL7G34350",1:10],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Root",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col=colnames(FPKM_mean)[1:10],angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,fontsize=12,
labels_row=c(expression(italic("FRO2")),expression(italic("COPT2")),expression(italic("bHLH39")),expression(italic("FRO4")),expression(italic("IRT1"))))$gtable
S<-pheatmap(FPKM_mean[rownames(FPKM_mean)=="AL1G10810"|rownames(FPKM_mean)=="AL5G25930"|rownames(FPKM_mean)=="AL5G38430"|rownames(FPKM_mean)=="AL6G35300"|rownames(FPKM_mean)=="AL7G34350",11:20],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Shoot",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col=colnames(FPKM_mean)[11:20],angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,fontsize=12,
labels_row=c(expression(italic("FRO2")),expression(italic("COPT2")),expression(italic("bHLH39")),expression(italic("FRO4")),expression(italic("IRT1"))))$gtable

g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Introgression_treatment_response.pdf",width=12,height=12,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()



IntroS1<-Kosi_S_Cu_Ath_desc_sorted_gff_lists[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists$Lyr_Gene),]
IntroS_KosiCtrl<-apply(IntroS1[,c(47:49)],1,mean)
IntroS_KosiCu<-apply(IntroS1[,c(44:46)],1,mean)
IntroS_KromCtrl<-apply(IntroS1[,c(41:43)],1,mean)
IntroS_KromCu<-apply(IntroS1[,c(38:40)],1,mean)
IntroS<-cbind(IntroS_KosiCtrl,IntroS_KosiCu,IntroS_KromCtrl,IntroS_KromCu)
rownames(IntroS)<-c("CCS","CV","CSD2","PKP4")
labels_row<-c(expression(italic(CCS)),expression(italic(CV)),expression(italic(CSD2)),expression(italic(PKP4)))
labels_col<-c("",expression(paste("Kosi ",bold(L),"Cu")),"","",expression(paste("Kosi ",bold(H),"Cu")),"","",expression(paste("Krom ",bold(L),"Cu")),"","",expression(paste("Krom ",bold(H),"Cu")),"")
IntroS<-IntroS[c(1,3,2,4),]
S<-pheatmap(IntroS, scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(Shoot),treeheight_row=0,treeheight_col=0,labels_col=labels_col,angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,labels_row=labels_row,fontsize=12)$gtable

require(gridExtra)
library(gtable)
library(grid)
g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Cu_reaction_stricter_beside_smaller.pdf",width=6,height=6,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()




