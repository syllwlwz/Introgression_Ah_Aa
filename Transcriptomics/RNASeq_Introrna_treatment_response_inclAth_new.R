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
#    R_1_1     R_2_1     R_3_1     R_4_1     R_5_1     R_6_1     R_7_1     R_8_1     R_9_1    R_10_1    R_32_1    R_33_2    R_34_1    R_35_1    R_44_1    R_45_1    R_46_1    R_49_1    R_50_1 
#1.0980109 1.0249802 1.0550967 1.1802253 1.1857587 1.0111051 0.8752349 0.9253697 0.9297485 0.9487734 1.0554862 0.9216191 1.0780057 0.8915333 1.1023014 0.9600079 0.9035737 1.0854878 1.1257056 
#   R_59_1    R_62_1    R_63_1    R_64_1    R_65_1    R_66_1    R_67_1    R_68_1    R_69_1    R_70_1    R_71_1 
#1.0764143 0.9279293 1.0037493 0.9769326 1.0447225 0.9645122 1.0064150 1.0132657 0.9761412 1.0071549 0.9206726 
dds_shoot.norm <-  estimateSizeFactors(dds_shoot)
sizeFactors(dds_shoot.norm)
#   S_12_1    S_14_1    S_16_1    S_18_1    S_20_1    S_22_1    S_24_1    S_26_1    S_28_1    S_30_1    S_36_1    S_38_1    S_40_1    S_42_1    S_47_1    S_51_1    S_53_1    S_55_1    S_57_1 
#1.0900708 0.9693277 1.0306529 1.0344110 1.0357065 0.8859311 0.9792407 1.1112269 1.0688144 0.9277719 0.9898966 0.9785856 0.9808552 1.1319682 0.9857119 1.0718547 0.9649453 1.0236305 1.1414299 
#   S_60_2    S_72_2    S_74_1    S_76_2    S_78_1    S_80_2    S_82_1    S_84_2    S_86_1    S_88_1    S_90_1 
#1.0214507 0.9540833 1.0121157 0.9854430 0.9596862 1.0257065 0.9059768 1.0104102 0.9473352 1.0692090 0.9896922 

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

# Performing estimation of dispersion parameter to account for different variances within populations for each gene
dds_root.disp <- estimateDispersions(dds_root.norm, fitType='local')
jpeg("Dispest_root_local.jpeg",width =8,height=8,units="in",res=1000)
plotDispEsts(dds_root.disp)
dev.off()

dds_root.disp <- estimateDispersions(dds_root.norm, fitType='parametric')
jpeg("Dispest_root_parametric.jpeg",width =8,height=8,units="in",res=1000)
plotDispEsts(dds_root.disp)
dev.off()

dds_root.disp <- estimateDispersions(dds_root.norm, fitType='mean')
jpeg("Dispest_root_mean.jpeg",width =8,height=8,units="in",res=1000)
plotDispEsts(dds_root.disp)
dev.off()

#go with default parametric, local automatically substituted if parametric does not converge,local gives more significant genes
#### Calculate Differential Expression
alpha <- 0.05 #significance level for adjusted p-value, default is 0.1

###############################################################################
#Differences between Populations and Species in treatment and treatment effect#
###############################################################################

design(dds_root)<-~ SpeciesPopcomb+Treatment+SpeciesPopcomb:Treatment
design(dds_shoot)<-~ SpeciesPopcomb+Treatment+SpeciesPopcomb:Treatment
dds_root_calc = DESeq(dds_root,test="Wald",fitType="parametric")
dds_shoot_calc = DESeq(dds_shoot,test="Wald",fitType="parametric")
resultsNames(dds_root_calc)
resultsNames(dds_shoot_calc)

#from the vignette:
#by adding genotype:condition, the main condition effect only represents the effect of condition for the reference level of genotype (I, or whichever level was defined by the user as the reference level). 
#The interaction terms genotypeII.conditionB and genotypeIII.conditionB give the difference between the condition effect for a given genotype and the condition effect for the reference genotype.

#from ?results:
## Example 3: two conditions, three genotypes

# ~~~ Using interaction terms ~~~

#dds <- makeExampleDESeqDataSet(n=100,m=18)
#dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
#design(dds) <- ~ genotype + condition + genotype:condition
#dds <- DESeq(dds)
#resultsNames(dds)
# the condition effect for genotype I (the main effect)
#should be Mias arenosa
IA_R_Treatment_arenosaMias<-results(dds_root_calc,contrast=c("Treatment","t","c"),alpha=0.05)
IA_R_Treatment_arenosaMias_sig<-as.data.frame(IA_R_Treatment_arenosaMias[which(IA_R_Treatment_arenosaMias$padj<=0.05&(IA_R_Treatment_arenosaMias$log2FoldChange>=1|IA_R_Treatment_arenosaMias$log2FoldChange<=(-1))),])
#97
IA_R_Treatment_arenosaMias_sig_0_5<-as.data.frame(IA_R_Treatment_arenosaMias[which(IA_R_Treatment_arenosaMias$padj<=0.05&(IA_R_Treatment_arenosaMias$log2FoldChange>=0.5|IA_R_Treatment_arenosaMias$log2FoldChange<=(-0.5))),])
#135

IA_S_Treatment_arenosaMias<-results(dds_shoot_calc,contrast=c("Treatment","t","c"),alpha=0.05)
IA_S_Treatment_arenosaMias_sig<-as.data.frame(IA_S_Treatment_arenosaMias[which(IA_S_Treatment_arenosaMias$padj<=0.05&(IA_S_Treatment_arenosaMias$log2FoldChange>=1|IA_S_Treatment_arenosaMias$log2FoldChange<=(-1))),])
#32
IA_S_Treatment_arenosaMias_sig_0_5<-as.data.frame(IA_S_Treatment_arenosaMias[which(IA_S_Treatment_arenosaMias$padj<=0.05&(IA_S_Treatment_arenosaMias$log2FoldChange>=0.5|IA_S_Treatment_arenosaMias$log2FoldChange<=(-0.5))),])
#35

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

#write.table(as.data.frame(FPKM_df),"FPKM_df_named.txt",sep="\t",row.names=FALSE)
FPKM_df<-read.table("FPKM_df.txt",sep="\t",header=T)

coldata_match<-coldata2[match(coldata2$Sample,colnames(FPKM_df)[3:62]),]
colnames(FPKM_df)[3:62]<-c(paste(coldata_match$Species,coldata_match$Population,coldata_match$Treatment,coldata_match$Tissue,coldata_match$Repeat,sep="_"))
colnames(FPKM_df)[63:122]<-paste("FPKM",colnames(FPKM_df)[3:62],sep="_")
countdata_FPKM<-data.frame(FPKM_df[,1:2],FPKM_df[,3:62][,order(colnames(FPKM_df[3:62]))],FPKM_df[,63:122][,order(colnames(FPKM_df[63:122]))])

IA_R_Treatment_arenosaMias_sig_lyr<-merge(IA_R_Treatment_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
#merge by rownames with thaliana gene description
Thal_description<-read.xlsx("E:/Lyrata_RBH/Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx")
All_lists<-read.xlsx("metal_homeostasis_20200707_v6sorted.xlsx",2)
All_lists<-All_lists[,-1]
All_lists$AGI_Number<-toupper(All_lists$AGI_Number)

IA_R_Treatment_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
#sort by fold change
IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted<-IA_R_Treatment_arenosaMias_sig_Ath_desc[order(-abs(IA_R_Treatment_arenosaMias_sig_Ath_desc$log2FoldChange)),]

#add gff description and info
gff<-read.table("D:/Lyrata/Alyrata_384_v2.1.gene.gff3",sep="\t",header=F)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
gff<-gff[gff$Type=="gene",]
ID<-substr(as.character(gff$Attributes),4,nchar(as.character(gff$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
gff$Lyr_Gene<-ID3[,1]

IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
#change order
IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted_gff<-IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]

#add lists
IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosaMias_inclAth.table",sep="\t",row.names=F)

IA_S_Treatment_arenosaMias_sig_lyr<-merge(IA_S_Treatment_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_Treatment_arenosaMias_sig_lyr)[1]<-"Alyr_ID"
IA_S_Treatment_arenosaMias_sig_Ath_desc<-merge(Thal_description,IA_S_Treatment_arenosaMias_sig_lyr,by="Alyr_ID",all.y=TRUE)
IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted<-IA_S_Treatment_arenosaMias_sig_Ath_desc[order(-abs(IA_S_Treatment_arenosaMias_sig_Ath_desc$log2FoldChange)),]
IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted_gff<-merge(IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted_gff<-IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted_gff_lists<-merge(IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_Treatment_arenosaMias_sig_Ath_desc_sorted_gff_lists,"IA_S_Treatment_arenosaMias_inclAth.table",sep="\t",row.names=F)

IA_R_Treatment_arenosaMias_sig_0_5_lyr<-merge(IA_R_Treatment_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosaMias_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosaMias_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosaMias_0_5_inclAth.table",sep="\t",row.names=F)

IA_S_Treatment_arenosaMias_sig_0_5_lyr<-merge(IA_S_Treatment_arenosaMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_Treatment_arenosaMias_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc<-merge(Thal_description,IA_S_Treatment_arenosaMias_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted<-IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc[order(-abs(IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc$log2FoldChange)),]
IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff<-merge(IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff<-IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff_lists,"IA_S_Treatment_arenosaMias_0_5_inclAth.table",sep="\t",row.names=F)

wb <- createWorkbook()
addWorksheet(wb, "Root")
addWorksheet(wb, "Shoot")

writeData(wb, "Root",IA_R_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Shoot",IA_S_Treatment_arenosaMias_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_MAa_treatment_lFC_0_5.xlsx", overwrite = TRUE)


# the condition effect for genotype III.
# this is the main effect *plus* the interaction term
# (the extra condition effect in genotype III compared to genotype I).
#results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))
 
IA_R_Treatment_arenosaZapa<-results(dds_root_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombarenosaZapa.Treatmentt")),alpha=0.05)
IA_R_Treatment_arenosaZapa_sig<-as.data.frame(IA_R_Treatment_arenosaZapa[which(IA_R_Treatment_arenosaZapa$padj<=0.05&(IA_R_Treatment_arenosaZapa$log2FoldChange>=1|IA_R_Treatment_arenosaZapa$log2FoldChange<=(-1))),])
#9
IA_R_Treatment_arenosaZapa_sig_0_5<-as.data.frame(IA_R_Treatment_arenosaZapa[which(IA_R_Treatment_arenosaZapa$padj<=0.05&(IA_R_Treatment_arenosaZapa$log2FoldChange>=0.5|IA_R_Treatment_arenosaZapa$log2FoldChange<=(-0.5))),])
#11

IA_S_Treatment_arenosaZapa<-results(dds_shoot_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombarenosaZapa.Treatmentt")),alpha=0.05)
IA_S_Treatment_arenosaZapa_sig<-as.data.frame(IA_S_Treatment_arenosaZapa[which(IA_S_Treatment_arenosaZapa$padj<=0.05&(IA_S_Treatment_arenosaZapa$log2FoldChange>=1|IA_S_Treatment_arenosaZapa$log2FoldChange<=(-1))),])
#1
IA_S_Treatment_arenosaZapa_sig_0_5<-as.data.frame(IA_S_Treatment_arenosaZapa[which(IA_S_Treatment_arenosaZapa$padj<=0.05&(IA_S_Treatment_arenosaZapa$log2FoldChange>=0.5|IA_S_Treatment_arenosaZapa$log2FoldChange<=(-0.5))),])
#1

IA_R_Treatment_arenosaZapa_sig_0_5_lyr<-merge(IA_R_Treatment_arenosaZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosaZapa_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosaZapa_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosaZapa_0_5_inclAth.table",sep="\t",row.names=F)

IA_S_Treatment_arenosaZapa_sig_0_5_lyr<-merge(IA_S_Treatment_arenosaZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_Treatment_arenosaZapa_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc<-merge(Thal_description,IA_S_Treatment_arenosaZapa_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted<-IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc[order(-abs(IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc$log2FoldChange)),]
IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff<-merge(IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff<-IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff_lists,"IA_S_Treatment_arenosaZapa_0_5_inclAth.table",sep="\t",row.names=F)

wb <- createWorkbook()
addWorksheet(wb, "Root")
addWorksheet(wb, "Shoot")

writeData(wb, "Root",IA_R_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Shoot",IA_S_Treatment_arenosaZapa_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_ZAa_treatment_lFC_0_5.xlsx", overwrite = TRUE)



IA_R_Treatment_halleriZapa<-results(dds_root_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombhalleriZapa.Treatmentt")),alpha=0.05)
IA_R_Treatment_halleriZapa_sig<-as.data.frame(IA_R_Treatment_halleriZapa[which(IA_R_Treatment_halleriZapa$padj<=0.05&(IA_R_Treatment_halleriZapa$log2FoldChange>=1|IA_R_Treatment_halleriZapa$log2FoldChange<=(-1))),])
#321
IA_R_Treatment_halleriZapa_sig_0_5<-as.data.frame(IA_R_Treatment_halleriZapa[which(IA_R_Treatment_halleriZapa$padj<=0.05&(IA_R_Treatment_halleriZapa$log2FoldChange>=0.5|IA_R_Treatment_halleriZapa$log2FoldChange<=(-0.5))),])
#660

IA_S_Treatment_halleriZapa<-results(dds_shoot_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombhalleriZapa.Treatmentt")),alpha=0.05)
IA_S_Treatment_halleriZapa_sig<-as.data.frame(IA_S_Treatment_halleriZapa[which(IA_S_Treatment_halleriZapa$padj<=0.05&(IA_S_Treatment_halleriZapa$log2FoldChange>=1|IA_S_Treatment_halleriZapa$log2FoldChange<=(-1))),])
#36
IA_S_Treatment_halleriZapa_sig_0_5<-as.data.frame(IA_S_Treatment_halleriZapa[which(IA_S_Treatment_halleriZapa$padj<=0.05&(IA_S_Treatment_halleriZapa$log2FoldChange>=0.5|IA_S_Treatment_halleriZapa$log2FoldChange<=(-0.5))),])
#39

IA_R_Treatment_halleriZapa_sig_0_5_lyr<-merge(IA_R_Treatment_halleriZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_halleriZapa_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_halleriZapa_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_halleriZapa_0_5_inclAth.table",sep="\t",row.names=F)

IA_S_Treatment_halleriZapa_sig_0_5_lyr<-merge(IA_S_Treatment_halleriZapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_Treatment_halleriZapa_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc<-merge(Thal_description,IA_S_Treatment_halleriZapa_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted<-IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc[order(-abs(IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc$log2FoldChange)),]
IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff<-merge(IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff<-IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff_lists,"IA_S_Treatment_halleriZapa_0_5_inclAth.table",sep="\t",row.names=F)

wb <- createWorkbook()
addWorksheet(wb, "Root")
addWorksheet(wb, "Shoot")

writeData(wb, "Root",IA_R_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Shoot",IA_S_Treatment_halleriZapa_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_ZAh_treatment_lFC_0_5.xlsx", overwrite = TRUE)



IA_R_Treatment_halleriMias<-results(dds_root_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombhalleriMias.Treatmentt")),alpha=0.05)
IA_R_Treatment_halleriMias_sig<-as.data.frame(IA_R_Treatment_halleriMias[which(IA_R_Treatment_halleriMias$padj<=0.05&(IA_R_Treatment_halleriMias$log2FoldChange>=1|IA_R_Treatment_halleriMias$log2FoldChange<=(-1))),])
#271
IA_R_Treatment_halleriMias_sig_0_5<-as.data.frame(IA_R_Treatment_halleriMias[which(IA_R_Treatment_halleriMias$padj<=0.05&(IA_R_Treatment_halleriMias$log2FoldChange>=0.5|IA_R_Treatment_halleriMias$log2FoldChange<=(-0.5))),])
#448

IA_S_Treatment_halleriMias<-results(dds_shoot_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombhalleriMias.Treatmentt")),alpha=0.05)
IA_S_Treatment_halleriMias_sig<-as.data.frame(IA_S_Treatment_halleriMias[which(IA_S_Treatment_halleriMias$padj<=0.05&(IA_S_Treatment_halleriMias$log2FoldChange>=1|IA_S_Treatment_halleriMias$log2FoldChange<=(-1))),])
#34
IA_S_Treatment_halleriMias_sig_0_5<-as.data.frame(IA_S_Treatment_halleriMias[which(IA_S_Treatment_halleriMias$padj<=0.05&(IA_S_Treatment_halleriMias$log2FoldChange>=0.5|IA_S_Treatment_halleriMias$log2FoldChange<=(-0.5))),])
#39

IA_R_Treatment_halleriMias_sig_0_5_lyr<-merge(IA_R_Treatment_halleriMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_halleriMias_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_halleriMias_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_halleriMias_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_halleriMias_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_halleriMias_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_halleriMias_0_5_inclAth.table",sep="\t",row.names=F)

IA_S_Treatment_halleriMias_sig_0_5_lyr<-merge(IA_S_Treatment_halleriMias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_Treatment_halleriMias_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_S_Treatment_halleriMias_sig_0_5_Ath_desc<-merge(Thal_description,IA_S_Treatment_halleriMias_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted<-IA_S_Treatment_halleriMias_sig_0_5_Ath_desc[order(-abs(IA_S_Treatment_halleriMias_sig_0_5_Ath_desc$log2FoldChange)),]
IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff<-merge(IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff<-IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff_lists,"IA_S_Treatment_halleriMias_0_5_inclAth.table",sep="\t",row.names=F)

wb <- createWorkbook()
addWorksheet(wb, "Root")
addWorksheet(wb, "Shoot")

writeData(wb, "Root",IA_R_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Shoot",IA_S_Treatment_halleriMias_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_MAh_treatment_lFC_0_5.xlsx", overwrite = TRUE)

IA_R_Treatment_thaliana<-results(dds_root_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombthalianaCol0.Treatmentt")),alpha=0.05)
IA_R_Treatment_thaliana_sig<-as.data.frame(IA_R_Treatment_thaliana[which(IA_R_Treatment_thaliana$padj<=0.05&(IA_R_Treatment_thaliana$log2FoldChange>=1|IA_R_Treatment_thaliana$log2FoldChange<=(-1))),])
#74
IA_R_Treatment_thaliana_sig_0_5<-as.data.frame(IA_R_Treatment_thaliana[which(IA_R_Treatment_thaliana$padj<=0.05&(IA_R_Treatment_thaliana$log2FoldChange>=0.5|IA_R_Treatment_thaliana$log2FoldChange<=(-0.5))),])
#94

IA_S_Treatment_thaliana<-results(dds_shoot_calc,contrast=list(c("Treatment_t_vs_c","SpeciesPopcombthalianaCol0.Treatmentt")),alpha=0.05)
IA_S_Treatment_thaliana_sig<-as.data.frame(IA_S_Treatment_thaliana[which(IA_S_Treatment_thaliana$padj<=0.05&(IA_S_Treatment_thaliana$log2FoldChange>=1|IA_S_Treatment_thaliana$log2FoldChange<=(-1))),])
#3
IA_S_Treatment_thaliana_sig_0_5<-as.data.frame(IA_S_Treatment_thaliana[which(IA_S_Treatment_thaliana$padj<=0.05&(IA_S_Treatment_thaliana$log2FoldChange>=0.5|IA_S_Treatment_thaliana$log2FoldChange<=(-0.5))),])
#3

IA_R_Treatment_thaliana_sig_0_5_lyr<-merge(IA_R_Treatment_thaliana_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_thaliana_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_thaliana_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_thaliana_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_thaliana_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_thaliana_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_thaliana_0_5_inclAth.table",sep="\t",row.names=F)

IA_S_Treatment_thaliana_sig_0_5_lyr<-merge(IA_S_Treatment_thaliana_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_Treatment_thaliana_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_S_Treatment_thaliana_sig_0_5_Ath_desc<-merge(Thal_description,IA_S_Treatment_thaliana_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted<-IA_S_Treatment_thaliana_sig_0_5_Ath_desc[order(-abs(IA_S_Treatment_thaliana_sig_0_5_Ath_desc$log2FoldChange)),]
IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff<-merge(IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff<-IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff_lists,"IA_S_Treatment_thaliana_0_5_inclAth.table",sep="\t",row.names=F)

wb <- createWorkbook()
addWorksheet(wb, "Root")
addWorksheet(wb, "Shoot")

writeData(wb, "Root",IA_R_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Shoot",IA_S_Treatment_thaliana_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_Ath_treatment_lFC_0_5.xlsx", overwrite = TRUE)




# the interaction term for condition effect in genotype III vs genotype I.
# this tests if the condition effect is different in III compared to I
#results(dds, name="genotypeIII.conditionB")

IA_R_Treatment_arenosa_Mias_vs_Zapa<-results(dds_root_calc,name="SpeciesPopcombarenosaZapa.Treatmentt",alpha=0.05)
IA_R_Treatment_arenosa_Mias_vs_Zapa_sig<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_Zapa[which(IA_R_Treatment_arenosa_Mias_vs_Zapa$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange>=1|IA_R_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange<=(-1))),])
#3
IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_Zapa[which(IA_R_Treatment_arenosa_Mias_vs_Zapa$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange>=0.5|IA_R_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange<=(-0.5))),])
#3

IA_S_Treatment_arenosa_Mias_vs_Zapa<-results(dds_shoot_calc,name="SpeciesPopcombarenosaZapa.Treatmentt",alpha=0.05)
IA_S_Treatment_arenosa_Mias_vs_Zapa_sig<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_Zapa[which(IA_S_Treatment_arenosa_Mias_vs_Zapa$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange>=1|IA_S_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange<=(-1))),])
#0
IA_S_Treatment_arenosa_Mias_vs_Zapa_sig_0_5<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_Zapa[which(IA_S_Treatment_arenosa_Mias_vs_Zapa$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange>=0.5|IA_S_Treatment_arenosa_Mias_vs_Zapa$log2FoldChange<=(-0.5))),])
#0

IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_lyr<-merge(IA_R_Treatment_arenosa_Mias_vs_Zapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosa_Mias_vs_Zapa_0_5_inclAth.table",sep="\t",row.names=F)

IA_R_Treatment_arenosa_Mias_vs_halleri_Mias<-results(dds_root_calc,name="SpeciesPopcombhalleriMias.Treatmentt",alpha=0.05)
IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias[which(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange>=1|IA_R_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange<=(-1))),])
#45
IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias[which(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange>=0.5|IA_R_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange<=(-0.5))),])
#52

IA_S_Treatment_arenosa_Mias_vs_halleri_Mias<-results(dds_shoot_calc,name="SpeciesPopcombhalleriMias.Treatmentt",alpha=0.05)
IA_S_Treatment_arenosa_Mias_vs_halleri_Mias_sig<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_halleri_Mias[which(IA_S_Treatment_arenosa_Mias_vs_halleri_Mias$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange>=1|IA_S_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange<=(-1))),])
#0
IA_S_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_halleri_Mias[which(IA_S_Treatment_arenosa_Mias_vs_halleri_Mias$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange>=0.5|IA_S_Treatment_arenosa_Mias_vs_halleri_Mias$log2FoldChange<=(-0.5))),])
#0

IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_lyr<-merge(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_0_5_inclAth.table",sep="\t",row.names=F)


IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa<-results(dds_root_calc,name="SpeciesPopcombhalleriZapa.Treatmentt",alpha=0.05)
IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa[which(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange>=1|IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange<=(-1))),])
#38
IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa[which(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange>=0.5|IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange<=(-0.5))),])
#49

IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa<-results(dds_shoot_calc,name="SpeciesPopcombhalleriZapa.Treatmentt",alpha=0.05)
IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa_sig<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa[which(IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange>=1|IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange<=(-1))),])
#0
IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa[which(IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange>=0.5|IA_S_Treatment_arenosa_Mias_vs_halleri_Zapa$log2FoldChange<=(-0.5))),])
#0

IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_lyr<-merge(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_0_5_inclAth.table",sep="\t",row.names=F)

IA_R_Treatment_arenosa_Mias_vs_thaliana<-results(dds_root_calc,name="SpeciesPopcombthalianaCol0.Treatmentt",alpha=0.05)
IA_R_Treatment_arenosa_Mias_vs_thaliana_sig<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_thaliana[which(IA_R_Treatment_arenosa_Mias_vs_thaliana$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange>=1|IA_R_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange<=(-1))),])
#17
IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5<-as.data.frame(IA_R_Treatment_arenosa_Mias_vs_thaliana[which(IA_R_Treatment_arenosa_Mias_vs_thaliana$padj<=0.05&(IA_R_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange>=0.5|IA_R_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange<=(-0.5))),])
#20

IA_S_Treatment_arenosa_Mias_vs_thaliana<-results(dds_shoot_calc,name="SpeciesPopcombthalianaCol0.Treatmentt",alpha=0.05)
IA_S_Treatment_arenosa_Mias_vs_thaliana_sig<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_thaliana[which(IA_S_Treatment_arenosa_Mias_vs_thaliana$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange>=1|IA_S_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange<=(-1))),])
#2
IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5<-as.data.frame(IA_S_Treatment_arenosa_Mias_vs_thaliana[which(IA_S_Treatment_arenosa_Mias_vs_thaliana$padj<=0.05&(IA_S_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange>=0.5|IA_S_Treatment_arenosa_Mias_vs_thaliana$log2FoldChange<=(-0.5))),])
#2

IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_lyr<-merge(IA_R_Treatment_arenosa_Mias_vs_thaliana_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosa_Mias_vs_thaliana_0_5_inclAth.table",sep="\t",row.names=F)

IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_lyr<-merge(IA_S_Treatment_arenosa_Mias_vs_thaliana_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc<-merge(Thal_description,IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted<-IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc[order(-abs(IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc$log2FoldChange)),]
IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff<-merge(IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff<-IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
write.table(IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff_lists,"IA_S_Treatment_arenosa_Mias_vs_thaliana_0_5_inclAth.table",sep="\t",row.names=F)


wb <- createWorkbook()
addWorksheet(wb, "Root_AaZ")
addWorksheet(wb, "Root_AhM")
addWorksheet(wb, "Root_AhZ")
addWorksheet(wb, "Root_Ath")
addWorksheet(wb, "Shoot_Ath")

writeData(wb, "Root_AaZ",IA_R_Treatment_arenosa_Mias_vs_Zapa_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Root_AhM",IA_R_Treatment_arenosa_Mias_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Root_AhZ",IA_R_Treatment_arenosa_Mias_vs_halleri_Zapa_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Root_Ath",IA_R_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)
writeData(wb, "Shoot_Ath",IA_S_Treatment_arenosa_Mias_vs_thaliana_sig_0_5_Ath_desc_sorted_gff_lists, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Deseq2_introgression_interactions_treatment_vs_AaM_lFC_0_5.xlsx", overwrite = TRUE)







# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))

# Note that a likelihood ratio could be used to test if there are any
# differences in the condition effect between the three genotypes.



##############################################################
#Zapa arenosa vs. Mias halleri#
#############################################################
dds_root2<-dds_root
dds_shoot2<-dds_shoot

dds_root2$SpeciesPopcomb<-relevel(dds_root2$SpeciesPopcomb,ref="arenosaZapa")
dds_shoot2$SpeciesPopcomb<-relevel(dds_shoot2$SpeciesPopcomb,ref="arenosaZapa")

design(dds_root2)<-~ SpeciesPopcomb+Treatment+SpeciesPopcomb:Treatment
design(dds_shoot2)<-~ SpeciesPopcomb+Treatment+SpeciesPopcomb:Treatment
dds_root2_calc = DESeq(dds_root2,test="Wald",fitType="parametric")
dds_shoot2_calc = DESeq(dds_shoot2,test="Wald",fitType="parametric")
resultsNames(dds_root2_calc)
resultsNames(dds_shoot2_calc)

#from the vignette:
#by adding genotype:condition, the main condition effect only represents the effect of condition for the reference level of genotype (I, or whichever level was defined by the user as the reference level). 
#The interaction terms genotypeII.conditionB and genotypeIII.conditionB give the difference between the condition effect for a given genotype and the condition effect for the reference genotype.

#from ?results:
## Example 3: two conditions, three genotypes

# ~~~ Using interaction terms ~~~

#dds <- makeExampleDESeqDataSet(n=100,m=18)
#dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
#design(dds) <- ~ genotype + condition + genotype:condition
#dds <- DESeq(dds)
#resultsNames(dds)
# the condition effect for genotype I (the main effect)
#should be Zapa arenosa
IA_R_Treatment_arenosaZapa<-results(dds_root2_calc,contrast=c("Treatment","t","c"),alpha=0.05)
IA_R_Treatment_arenosaZapa_sig<-as.data.frame(IA_R_Treatment_arenosaZapa[which(IA_R_Treatment_arenosaZapa$padj<=0.05&(IA_R_Treatment_arenosaZapa$log2FoldChange>=1|IA_R_Treatment_arenosaMias$log2FoldChange<=(-1))),])
#97
IA_R_Treatment_arenosaMias_sig_0_5<-as.data.frame(IA_R_Treatment_arenosaMias[which(IA_R_Treatment_arenosaMias$padj<=0.05&(IA_R_Treatment_arenosaMias$log2FoldChange>=0.5|IA_R_Treatment_arenosaMias$log2FoldChange<=(-0.5))),])
#135

IA_S_Treatment_arenosaMias<-results(dds_shoot_calc,contrast=c("Treatment","t","c"),alpha=0.05)
IA_S_Treatment_arenosaMias_sig<-as.data.frame(IA_S_Treatment_arenosaMias[which(IA_S_Treatment_arenosaMias$padj<=0.05&(IA_S_Treatment_arenosaMias$log2FoldChange>=1|IA_S_Treatment_arenosaMias$log2FoldChange<=(-1))),])
#32
IA_S_Treatment_arenosaMias_sig_0_5<-as.data.frame(IA_S_Treatment_arenosaMias[which(IA_S_Treatment_arenosaMias$padj<=0.05&(IA_S_Treatment_arenosaMias$log2FoldChange>=0.5|IA_S_Treatment_arenosaMias$log2FoldChange<=(-0.5))),])
#35





IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias<-results(dds_root_calc,name="SpeciesPopcombhalleriMias.Treatmentt",alpha=0.05)
IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig<-as.data.frame(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias[which(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias$padj<=0.05&(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange>=1|IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange<=(-1))),])
#45
IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5<-as.data.frame(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias[which(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias$padj<=0.05&(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange>=0.5|IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange<=(-0.5))),])
#52

IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias<-results(dds_shoot_calc,name="SpeciesPopcombhalleriMias.Treatmentt",alpha=0.05)
IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias_sig<-as.data.frame(IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias[which(IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias$padj<=0.05&(IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange>=1|IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange<=(-1))),])
#0
IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5<-as.data.frame(IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias[which(IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias$padj<=0.05&(IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange>=0.5|IA_S_Treatment_arenosa_Zapa_vs_halleri_Mias$log2FoldChange<=(-0.5))),])
#0

IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_lyr<-merge(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig,countdata_FPKM,by.x="row.names",by.y="Alyr_ID",all.x=TRUE)
colnames(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_lyr)[1]<-"Alyr_ID"
IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc<-merge(Thal_description,IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_lyr,by="Alyr_ID",all.y=TRUE)
IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted<-IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc[order(-abs(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc$log2FoldChange)),]
IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff<-merge(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted,gff,by.x="Alyr_ID",by.y="Lyr_Gene",all.x=T)
colnames(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff)[141:149]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff<-IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff[,c(1,141:149,2:140)]
IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff_lists<-merge(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
#write.table(IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_sig_0_5_Ath_desc_sorted_gff_lists,"IA_R_Treatment_arenosa_Zapa_vs_halleri_Mias_0_5_inclAth.table",sep="\t",row.names=F)





















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

#AL7G22090 HMA2
#AL1G64230 HAC04
#AL3G48300

require(gridExtra)
library(gtable)
library(grid)

R<-pheatmap(FPKM_mean[rownames(FPKM_mean)=="AL7G22090"|rownames(FPKM_mean)=="AL1G64230"|rownames(FPKM_mean)=="AL3G48300",1:10],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Root",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col=colnames(FPKM_mean)[1:10],angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,fontsize=12,
labels_row=c(expression(italic("HAC04")),expression(italic("unknown")),expression(italic("HMA2"))))$gtable

grid.newpage()
pdf("Heatmap_Introgression_treatment_interaction_Mias_Zapa_arenosa.pdf",width=12,height=12,paper="special")
grid.arrange(R)
dev.off()



#AL5G25930 COPT2
#AL4G15434 GLR2.3_1
#AL4G15440 GLR2.3_2
#AL6G24350 ZIF1
#AL6G35310 FRO5


R<-pheatmap(FPKM_mean[rownames(FPKM_mean)=="AL5G25930"|rownames(FPKM_mean)=="AL4G15434"|rownames(FPKM_mean)=="AL4G15440"|rownames(FPKM_mean)=="AL6G24350"|rownames(FPKM_mean)=="AL6G35310",1:10],scale="row",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F,cluster_cols=F,main="Root",treeheight_row =0,treeheight_col=0,show_colnames=T,labels_col=colnames(FPKM_mean)[1:10],angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,fontsize=12,
labels_row=c(expression(italic("GLR2.3_1")),expression(italic("GLR2.3_2")),expression(italic("COPT2")),expression(italic("ZIF1")),expression(italic("FRO5"))))$gtable

grid.newpage()
pdf("Heatmap_HM_Introgression_treatment_interaction_Mias_arenosa_halleri.pdf",width=12,height=12,paper="special")
grid.arrange(R)
dev.off()












