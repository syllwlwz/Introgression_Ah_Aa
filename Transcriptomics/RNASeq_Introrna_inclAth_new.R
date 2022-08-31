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

#correct HMA4
countdata2[rownames(countdata2)=="AL3G52820",1:60]<-countdata2[rownames(countdata2)=="AL3G52820",1:60]+countdata2[rownames(countdata2)=="AL1G54170",1:60]
countdata2<-countdata2[!rownames(countdata2)=="AL1G54170",]

#31072 instead of 31073 genes

countdata[rownames(countdata)=="AL3G52820",1:60]<-countdata[rownames(countdata)=="AL3G52820",1:60]+countdata[rownames(countdata)=="AL1G54170",1:60]
countdata<-countdata[!rownames(countdata)=="AL1G54170",]



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
dds2 <- DESeqDataSetFromMatrix(countData=countdata2,colData=coldata2,design=~Tissue+Treatment+SpeciesPopcomb)

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

colData(dds_root)<-colData(dds_root)[,c(1:7)]
colData(dds_shoot)<-colData(dds_shoot)[,c(1:7)]

#Note: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the
#formula and make sure the control level is the first level.

#### Normalization
# Normalizing for different numbers of aligned reads per library  
dds_root.norm <-  estimateSizeFactors(dds_root)
sizeFactors(dds_root.norm)
#   R_1_1     R_2_1     R_3_1     R_4_1     R_5_1     R_6_1     R_7_1     R_8_1 
#1.0979910 1.0249857 1.0551706 1.1802265 1.1858042 1.0110724 0.8751955 0.9253561 
#    R_9_1    R_10_1    R_32_1    R_33_2    R_34_1    R_35_1    R_44_1    R_45_1 
#0.9297359 0.9487662 1.0555204 0.9216135 1.0779914 0.8915570 1.1022784 0.9600023 
#   R_46_1    R_49_1    R_50_1    R_59_1    R_62_1    R_63_1    R_64_1    R_65_1 
#0.9035927 1.0854841 1.1257247 1.0763884 0.9279306 1.0037240 0.9769464 1.0447017 
#   R_66_1    R_67_1    R_68_1    R_69_1    R_70_1    R_71_1 
#0.9644947 1.0063855 1.0132534 0.9761422 1.0072074 0.9206471
dds_shoot.norm <-  estimateSizeFactors(dds_shoot)
sizeFactors(dds_shoot.norm)
#   S_12_1    S_14_1    S_16_1    S_18_1    S_20_1    S_22_1    S_24_1    S_26_1 
#1.0900377 0.9692911 1.0306531 1.0344272 1.0356918 0.8859463 0.9792458 1.1111969 
#   S_28_1    S_30_1    S_36_1    S_38_1    S_40_1    S_42_1    S_47_1    S_51_1 
#1.0688137 0.9277488 0.9899324 0.9785351 0.9807994 1.1319819 0.9857152 1.0718508 
#   S_53_1    S_55_1    S_57_1    S_60_2    S_72_2    S_74_1    S_76_2    S_78_1 
#0.9649672 1.0236228 1.1414229 1.0214496 0.9541365 1.0119177 0.9854292 0.9596967 
#   S_80_2    S_82_1    S_84_2    S_86_1    S_88_1    S_90_1 
#1.0256933 0.9059744 1.0104007 0.9473343 1.0692150 0.9897604

require(affy)
# Checking the normalization
jpeg("Distributions_raw_counts_and_normalized_counts_root_inclAth.jpeg",width =8.447917,height=6.072917,units="in",res=1000)
par(mfrow=c(2,2),mar=c(2,7,1,0.5),cex.lab=0.7)#set parameters for the plotting window
epsilon <- 1 # pseudo-count to avoid problems with log(0)
boxplot(log2(counts(dds_root.norm)+epsilon),col=c("red","red","red","grey","grey","grey"), cex.axis=0.7, 
        las=1, xlab="log2(counts+1)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds_root.norm, normalized=TRUE)+epsilon), col=c("red","red","red","grey","grey","grey"),cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds_root.norm)+epsilon), 
            xlab="log2(counts+1)", col=c("red","red","red","grey","grey","grey"),cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds_root.norm, normalized=TRUE)+epsilon), 
            xlab="log2(normalized counts)",col=c("red","red","red","grey","grey","grey"), cex.lab=0.7, panel.first=grid()) 
dev.off()
jpeg("Distributions_raw_counts_and_normalized_counts_shoot_inclAth.jpeg",width =8.447917,height=6.072917,units="in",res=1000)
par(mfrow=c(2,2),mar=c(2,7,1,0.5),cex.lab=0.7)#set parameters for the plotting window
epsilon <- 1 # pseudo-count to avoid problems with log(0)
boxplot(log2(counts(dds_shoot.norm)+epsilon),col=c("red","red","red","grey","grey","grey"), cex.axis=0.7, 
        las=1, xlab="log2(counts+1)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds_shoot.norm, normalized=TRUE)+epsilon), col=c("red","red","red","grey","grey","grey"),cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds_shoot.norm)+epsilon), 
            xlab="log2(counts+1)", col=c("red","red","red","grey","grey","grey"),cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds_shoot.norm, normalized=TRUE)+epsilon), 
            xlab="log2(normalized counts)",col=c("red","red","red","grey","grey","grey"), cex.lab=0.7, panel.first=grid()) 
dev.off()

# Restore default parameters
par(mfrow=c(1,1), cex.lab=1,mar=c(5.1, 4.1, 4.1, 2.1))

#go with default parametric, local automatically substituted if parametric does not converge,local gives more significant genes
#### Calculate Differential Expression
alpha <- 0.05 #significance level for adjusted p-value, default is 0.1

######################################################################
#Differences between Populations and Species independent of treatment#
######################################################################

design(dds_root)<-~ SpeciesPopcomb+Treatment+SpeciesPopcomb:Treatment
design(dds_shoot)<-~ SpeciesPopcomb+Treatment+SpeciesPopcomb:Treatment
dds_root_calc = DESeq(dds_root,test="Wald",fitType="parametric")
dds_shoot_calc = DESeq(dds_shoot,test="Wald",fitType="parametric")
resultsNames(dds_root_calc)
resultsNames(dds_shoot_calc)

IA_R_arenosaZapa_vs_arenosaMias<-results(dds_root_calc,contrast=c("SpeciesPopcomb","arenosaZapa","arenosaMias"),alpha=0.05)
IA_R_halleriMias_vs_arenosaMias<-results(dds_root_calc,contrast=c("SpeciesPopcomb","halleriMias","arenosaMias"),alpha=0.05)
IA_S_arenosaZapa_vs_arenosaMias<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","arenosaZapa","arenosaMias"),alpha=0.05)
IA_S_halleriMias_vs_arenosaMias<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","halleriMias","arenosaMias"),alpha=0.05)
IA_R_halleriZapa_vs_arenosaMias<-results(dds_root_calc,contrast=c("SpeciesPopcomb","halleriZapa","arenosaMias"),alpha=0.05)
IA_S_halleriZapa_vs_arenosaMias<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","halleriZapa","arenosaMias"),alpha=0.05)
IA_R_thalianaCol0_vs_arenosaMias<-results(dds_root_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaMias"),alpha=0.05)
IA_S_thalianaCol0_vs_arenosaMias<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaMias"),alpha=0.05)
IA_R_thalianaCol0_vs_arenosaZapa<-results(dds_root_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaZapa"),alpha=0.05)
IA_S_thalianaCol0_vs_arenosaZapa<-results(dds_shoot_calc,contrast=c("SpeciesPopcomb","thalianaCol0","arenosaZapa"),alpha=0.05)

IA_R_arenosaZapa_vs_arenosaMias_sig<-as.data.frame(IA_R_arenosaZapa_vs_arenosaMias[which(IA_R_arenosaZapa_vs_arenosaMias$padj<=0.05&(IA_R_arenosaZapa_vs_arenosaMias$log2FoldChange>=1|IA_R_arenosaZapa_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_R_halleriMias_vs_arenosaMias_sig<-as.data.frame(IA_R_halleriMias_vs_arenosaMias[which(IA_R_halleriMias_vs_arenosaMias$padj<=0.05&(IA_R_halleriMias_vs_arenosaMias$log2FoldChange>=1|IA_R_halleriMias_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_S_arenosaZapa_vs_arenosaMias_sig<-as.data.frame(IA_S_arenosaZapa_vs_arenosaMias[which(IA_S_arenosaZapa_vs_arenosaMias$padj<=0.05&(IA_S_arenosaZapa_vs_arenosaMias$log2FoldChange>=1|IA_S_arenosaZapa_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_S_halleriMias_vs_arenosaMias_sig<-as.data.frame(IA_S_halleriMias_vs_arenosaMias[which(IA_S_halleriMias_vs_arenosaMias$padj<=0.05&(IA_S_halleriMias_vs_arenosaMias$log2FoldChange>=1|IA_S_halleriMias_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_R_halleriZapa_vs_arenosaMias_sig<-as.data.frame(IA_R_halleriZapa_vs_arenosaMias[which(IA_R_halleriZapa_vs_arenosaMias$padj<=0.05&(IA_R_halleriZapa_vs_arenosaMias$log2FoldChange>=1|IA_R_halleriZapa_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_S_halleriZapa_vs_arenosaMias_sig<-as.data.frame(IA_S_halleriZapa_vs_arenosaMias[which(IA_S_halleriZapa_vs_arenosaMias$padj<=0.05&(IA_S_halleriZapa_vs_arenosaMias$log2FoldChange>=1|IA_S_halleriZapa_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_R_thalianaCol0_vs_arenosaMias_sig<-as.data.frame(IA_R_thalianaCol0_vs_arenosaMias[which(IA_R_thalianaCol0_vs_arenosaMias$padj<=0.05&(IA_R_thalianaCol0_vs_arenosaMias$log2FoldChange>=1|IA_R_thalianaCol0_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_S_thalianaCol0_vs_arenosaMias_sig<-as.data.frame(IA_S_thalianaCol0_vs_arenosaMias[which(IA_S_thalianaCol0_vs_arenosaMias$padj<=0.05&(IA_S_thalianaCol0_vs_arenosaMias$log2FoldChange>=1|IA_S_thalianaCol0_vs_arenosaMias$log2FoldChange<=(-1))),])
IA_R_thalianaCol0_vs_arenosaZapa_sig<-as.data.frame(IA_R_thalianaCol0_vs_arenosaZapa[which(IA_R_thalianaCol0_vs_arenosaZapa$padj<=0.05&(IA_R_thalianaCol0_vs_arenosaZapa$log2FoldChange>=1|IA_R_thalianaCol0_vs_arenosaZapa$log2FoldChange<=(-1))),])
IA_S_thalianaCol0_vs_arenosaZapa_sig<-as.data.frame(IA_S_thalianaCol0_vs_arenosaZapa[which(IA_S_thalianaCol0_vs_arenosaZapa$padj<=0.05&(IA_S_thalianaCol0_vs_arenosaZapa$log2FoldChange>=1|IA_S_thalianaCol0_vs_arenosaZapa$log2FoldChange<=(-1))),])

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

write.table(as.data.frame(FPKM_df),"FPKM_df_named.txt",sep="\t",row.names=FALSE)
#FPKM_df<-read.table("FPKM_df.txt",sep="\t",header=T)

#separate root and shoot
dds_root_norm <- estimateSizeFactors(dds_root)
normalized_counts_R <-data.frame(counts(dds_root_norm, normalized=TRUE))
colnames(normalized_counts_R)<-dds_root_norm$Sample.Name

dds_shoot_norm <- estimateSizeFactors(dds_shoot)
normalized_counts_S <-data.frame(counts(dds_shoot_norm, normalized=TRUE))
colnames(normalized_counts_S)<-dds_shoot_norm$Sample.Name

countdata_df<-cbind(normalized_counts_R,normalized_counts_S)
Extracted_coldata_R<-as.data.frame(colData(dds_root))
Extracted_coldata_S<-as.data.frame(colData(dds_shoot))

colnames(countdata_df)<-c(paste(Extracted_coldata_R$Species,Extracted_coldata_R$Population,Extracted_coldata_R$Treatment,Extracted_coldata_R$Tissue,Extracted_coldata_R$Repeat,sep="_"),paste(Extracted_coldata_S$Species,Extracted_coldata_S$Population,Extracted_coldata_S$Treatment,Extracted_coldata_S$Tissue,Extracted_coldata_S$Repeat,sep="_"))
countdata_rlog<-countdata_df[,order(colnames(countdata_df))]

countdata_rlog2<-merge(CDS_df[,c(1,3)],countdata_rlog,by.x="seq_name",by.y="row.names",all.y=T)

countdata_rlog2[,3:50]<-countdata_rlog2[,3:50]*1000/countdata_rlog2$Tlen
countdata_rlog_df<-countdata_rlog2[,-2]
countdata_rlog_df_Tair<-merge(countdata_rlog_df,Thal_description,by.x="seq_name",by.y="Alyr_ID",all.x=TRUE)

write.table(countdata_rlog_df_Tair,"Countdata_medofratios.table",sep="\t",row.names=F)


#count0 = cbind(genes2, datafrFC) #merge with gene name
#countdata <- count0[,-1]
#rownames(countdata) <- count0[,1] 
#colnames(countdata) <- gsub(".counts.txt","",gsub("./Counts_FC/","",filelistFC)) # add column names
#mode(countdata)
#mode(countdata) = "numeric" #convert the class of a matrix from character to numeric
#class(countdata)

#statsPerSample <- data.frame(t(apply(countdata, 2, summary)))
#head(statsPerSample)

coldata_match<-coldata2[match(coldata2$Sample,colnames(FPKM_df)[3:62]),]
colnames(FPKM_df)[3:62]<-c(paste(coldata_match$Species,coldata_match$Population,coldata_match$Treatment,coldata_match$Tissue,coldata_match$Repeat,sep="_"))
colnames(FPKM_df)[63:122]<-paste("FPKM",colnames(FPKM_df)[3:62],sep="_")

countdata_FPKM<-data.frame(FPKM_df[,1:2],FPKM_df[,3:62][,order(colnames(FPKM_df[3:62]))],FPKM_df[,63:122][,order(colnames(FPKM_df[63:122]))])
write.table(as.data.frame(countdata_FPKM),"FPKM_df_named.txt",sep="\t",row.names=FALSE)


IA_R_arenosaZapa_vs_arenosaMias_sig_lyr<-merge(IA_R_arenosaZapa_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_arenosaZapa_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_R_arenosaZapa_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted<-IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

#add gff description and info
gff<-read.table("D:/Lyrata/Alyrata_384_v2.1.gene.gff3",sep="\t",header=F)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
gff<-gff[gff$Type=="gene",]
ID<-substr(as.character(gff$Attributes),4,nchar(as.character(gff$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
gff$Lyr_Gene<-ID3[,1]

IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_R_arenosaZapa_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)

IA_S_arenosaZapa_vs_arenosaMias_sig_lyr<-merge(IA_S_arenosaZapa_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_arenosaZapa_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_S_arenosaZapa_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted<-IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_S_arenosaZapa_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)

IA_R_halleriMias_vs_arenosaMias_sig_lyr<-merge(IA_R_halleriMias_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_halleriMias_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_R_halleriMias_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted<-IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_R_halleriMias_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)


IA_S_halleriMias_vs_arenosaMias_sig_lyr<-merge(IA_S_halleriMias_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_halleriMias_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_S_halleriMias_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted<-IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_S_halleriMias_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)


IA_R_halleriZapa_vs_arenosaMias_sig_lyr<-merge(IA_R_halleriZapa_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_halleriZapa_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_R_halleriZapa_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted<-IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_R_halleriZapa_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)


IA_S_halleriZapa_vs_arenosaMias_sig_lyr<-merge(IA_S_halleriZapa_vs_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_halleriZapa_vs_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_S_halleriZapa_vs_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted<-IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc[order(-abs(IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc$log2FoldChange)),]

IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff<-IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_S_halleriZapa_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)


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
write.table(IA_R_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_R_thalianaCol0_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)


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
write.table(IA_S_thalianaCol0_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_S_thalianaCol0_vs_arenosaMias_inclAth.table",sep="\t",row.names=F)


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
write.table(IA_R_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff_lists,"IA_R_thalianaCol0_vs_arenosaZapa_inclAth.table",sep="\t",row.names=F)


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
write.table(IA_S_thalianaCol0_vs_arenosaZapa_sig_Ath_desc_sorted_gff_lists,"IA_S_thalianaCol0_vs_arenosaZapa_inclAth.table",sep="\t",row.names=F)


###################################################################################
#sign. diff in Mias arenosa vs. Zapa arenosa but not Mias halleri vs. Mias arenosa#
###################################################################################
Cand_R<-IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists[!IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID%in%IA_R_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID,]
#1976 of 4375
Cand_S<-IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists[!IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID%in%IA_S_halleriMias_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID,]
#1306 of 2901

wb <- createWorkbook()
addWorksheet(wb, "Root")
addWorksheet(wb, "Shoot")

writeData(wb, "Root",Cand_R, startRow = 1, startCol = 1)
writeData(wb, "Shoot",Cand_S, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_inclAth.xlsx", overwrite = TRUE)


###################################################################################
#sign. diff in Mias arenosa vs. Zapa arenosa but not Zapa halleri vs. Mias arenosa#
###################################################################################

Cand_R2<-IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists[!IA_R_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID%in%IA_R_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID,]
#1816 of 4375
Cand_S2<-IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists[!IA_S_arenosaZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID%in%IA_S_halleriZapa_vs_arenosaMias_sig_Ath_desc_sorted_gff_lists$Alyr_ID,]
#1283 of 2901

wb <- createWorkbook()
addWorksheet(wb, "Root")
addWorksheet(wb, "Shoot")

writeData(wb, "Root",Cand_R2, startRow = 1, startCol = 1)
writeData(wb, "Shoot",Cand_S2, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_MAa_vs_ZAa_wo_MAa_vs_ZAh_inclAth.xlsx", overwrite = TRUE)


















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

pdf("PCA_all_introgression_new.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA(Log, intgroup=c("Population", "Tissue","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=as.factor(pcaData$Tissue),pch=c(16,17,15)[as.factor(pcaData$Species)],xlab="",ylab="",main="Total",cex=1.5,cex.lab=1.2,axes=F)
legend("topright",inset=c(-0.22,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,pch=c(16,17,15),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Root","Shoot"),fill=c("black","red"),title="Tissue",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))
pcaData <- plotPCA(LogR, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=c("green","red","black")[as.factor(pcaData$Population)],pch=c(16,17,15)[as.factor(pcaData$Species)],xlab="",ylab="",main="Root",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA(LogS, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31073)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=c("green","red","black")[as.factor(pcaData$Population)],pch=c(16,17,15)[as.factor(pcaData$Species)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.5,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,pch=c(16,17,15),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.35,0.3),legend=c("Col0","Mias","Zapa"),fill=c("green","red","black"),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)

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



pdf("PCA_all_introgression_PC3_4_new.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC3_4(Log, intgroup=c("Population", "Tissue","Species","Treatment"), returnData=TRUE,ntop=31072)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC4~pcaData$PC3,col=as.factor(pcaData$Tissue),pch=c(16,17,15)[as.factor(pcaData$Species)],xlab="",ylab="",main="Total",cex=1.5,cex.lab=1.2,axes=F)
legend("topright",inset=c(-0.22,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,pch=c(16,17,15),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Root","Shoot"),fill=c("black","red"),title="Tissue",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC3 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC4 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))
pcaData <- plotPCA.PC3_4(LogR, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31072)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC4~pcaData$PC3,col=c("green","red","black")[as.factor(pcaData$Population)],pch=c(16,17,15)[as.factor(pcaData$Species)],xlab="",ylab="",main="Root",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC3 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC4 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA.PC3_4(LogS, intgroup=c("Population","Species","Treatment"), returnData=TRUE,ntop=31072)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC4~pcaData$PC3,col=c("green","red","black")[as.factor(pcaData$Population)],pch=c(16,17,15)[as.factor(pcaData$Species)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC3 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC4 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.5,0),legend=c("A. arenosa","A. halleri","A. thaliana"),text.font=3,pch=c(16,17,15),title="Species",bty="n",cex=1.2)
legend("topright",inset=c(-0.35,0.3),legend=c("Col0","Mias","Zapa"),fill=c("green","red","black"),title="Population",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)

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

#AL3G52820 HMA4
#AL4G46270 MTP1
#AL7G22080 HMA3
#AL4G24850 ZIP6
#AL4G13010 NRAMP3
#AL1G64050 ZIP11
Cand<-c("AL3G52820","AL4G46270","AL7G22080","AL4G24850","AL4G13010","AL1G64050")

LogR_assay<-assay(LogR)
colnames(LogR_assay)<-paste(LogR@colData@ listData$Species,LogR@colData@ listData$Population,LogR@colData@ listData$Treatment,sep="-" )
LogR_mean<-t(apply(LogR_assay,1,function(x) tapply(x,colnames(LogR_assay),mean)))

Cand_LogR<-LogR_mean[rownames(LogR_mean)%in%Cand,1:10]

LogS_assay<-assay(LogS)
colnames(LogS_assay)<-paste(LogS@colData@ listData$Species,LogS@colData@ listData$Population,LogS@colData@ listData$Treatment,sep="-" )
LogS_mean<-t(apply(LogS_assay,1,function(x) tapply(x,colnames(LogS_assay),mean)))

Cand_LogS<-LogS_mean[rownames(LogS_mean)%in%Cand,1:10]

R<-pheatmap(Cand_LogR[c(2,5,6,4,3,1),],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Root",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col="",angle_col=45,fontsize_col=10,cellwidth=10,cellheight=10,legend=F,fontsize=8,fontsize_row=8,
labels_row="")$gtable
S<-pheatmap(Cand_LogS[c(2,5,6,4,3,1),],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Shoot",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col="",angle_col=45,fontsize_col=10,cellwidth=10,cellheight=10,legend=F,fontsize=8,fontsize_row=8,
labels_row=c(expression(italic("HMA4")),expression(italic("MTP1")),expression(italic("HMA3")),expression(italic("ZIP6")),expression(italic("NRAMP3")),expression(italic("ZIP11"))))$gtable


require(gridExtra)
library(gtable)
library(grid)

g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Introgression_HM_core_cand_rlog.pdf",width=6,height=6,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()


#AL7G34360 IRT2
#AL8G14120 ZIP8
#AL6G35310 FRO5
#AL6G10450 FER1
#AL3G19430 FRD3

Cand<-c("AL7G34360","AL8G14120","AL6G35310","AL6G10450","AL3G19430")

Cand_LogR<-LogR_mean[rownames(LogR_mean)%in%Cand,1:10]

Cand_LogS<-LogS_mean[rownames(LogS_mean)%in%Cand,1:10]

R<-pheatmap(Cand_LogR[c(4,5,3,2,1),],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Root",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col="",angle_col=45,fontsize_col=10,cellwidth=10,cellheight=10,legend=F,fontsize=8,fontsize_row=8,
labels_row="")$gtable
S<-pheatmap(Cand_LogS[c(4,5,3,2,1),],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Shoot",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col="",angle_col=45,fontsize_col=10,cellwidth=10,cellheight=10,legend=F,fontsize=8,fontsize_row=8,
labels_row=c(expression(italic("IRT2")),expression(italic("ZIP8")),expression(italic("FRO5")),expression(italic("FER1")),expression(italic("FRD3"))))$gtable


require(gridExtra)
library(gtable)
library(grid)

g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Introgression_HM_noncore_cand_rlog.pdf",width=6,height=6,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()


#AL1G64230 HAC04
#AL7G22090 HMA2
#AL3G48300


Cand<-c("AL1G64230","AL7G22090","AL3G48300")

Cand_LogR<-LogR_mean[rownames(LogR_mean)%in%Cand,1:10]

Cand_LogS<-LogS_mean[rownames(LogS_mean)%in%Cand,1:10]

R<-pheatmap(Cand_LogR[c(3,1,2),],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Root",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col="",angle_col=45,fontsize_col=10,cellwidth=10,cellheight=10,legend=F,fontsize=8,fontsize_row=8,
labels_row="")$gtable
S<-pheatmap(Cand_LogS[c(3,1,2),],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Shoot",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col="",angle_col=45,fontsize_col=10,cellwidth=10,cellheight=10,legend=F,fontsize=8,fontsize_row=8,
labels_row=c(expression(italic("HMA2")),expression(italic("HAC04")),expression(italic("unknown"))))$gtable


require(gridExtra)
library(gtable)
library(grid)

g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Introgression_HM_interaction_treatment_cand_rlog.pdf",width=6,height=6,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()






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
#AL7G22080 HMA3
#AL4G24850 ZIP6
#AL4G13010 NRAMP3
#AL1G64050 ZIP11



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




