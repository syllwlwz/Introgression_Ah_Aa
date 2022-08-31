---
title: "Twisst"
author: "Lara Syllwasschy"
date: "12 October 2021"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
#modified from
#---
#title: "Twisst_analysis"
#author: "Patrick Monnahan"
#date: "12/29/2018"
#output: html_document
#---

##Load data and packages
```{r setup, include=FALSE}
setwd("E:/Introgression/Twisst/MiashalleriKowar_w100")
knitr::opts_chunk$set(echo = TRUE)
chromNames = paste0("scaffold_",1:8)
weights_files4 = paste0("MiashalleriKowaarenosa.weights.w100.m1.",chromNames,".csv") 
window_data_files = paste0("MiashalleriKowaarenosalyr.Vcf.scf",chromNames,".w100m1.LDphase.phyml_bionj.data.tsv")
#make list of weight and window data for each chromsome, throwing out windows with NA values
weights_by_chrom4 = list()
window_data_by_chrom4 = list()
for (i in 1:length(weights_files4)){
    weights <- read.table(weights_files4[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom4[[i]] <- weights[good_rows,]
    window_data_by_chrom4[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }
weights_all4=data.frame()
for(i in 1:8){
weights_all4=rbind(weights_all4,weights_by_chrom4[[i]])}
mean_weights4 <- apply(weights_all4, 2, mean)
twisst4=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom4[[i]],weights_by_chrom4[[i]])
    twisst4=rbind(twisst4,df)
}

require(ape)
pdf("Twisst_MiashalleriKowaar_weightings_w100_simple.pdf",width=8,height=8,paper="special",pointsize=20)
layout(matrix(c(1,2,3,4,5,5,5,5,5,5,5,5), nrow=3, byrow=TRUE),widths=c(0.5,1,1,1,1))
#topo1 (A_ly,(Mias_arenosa,(halleri_halleri,Kowa_arenosa)));
#topo2 (A_ly,(halleri_halleri,(Mias_arenosa,Kowa_arenosa)));
#topo3 (A_ly,(Kowa_arenosa,(Mias_arenosa,halleri)));

topo1 <-ape::read.tree(text='(A_ly,(Mias_arenosa,(halleri,Kowa_arenosa)));')
topo2<- ape::read.tree(text='(A_ly,(halleri,(Mias_arenosa,Kowa_arenosa)));')
topo3<- ape::read.tree(text='(A_ly,(Kowa_arenosa,(Mias_arenosa,halleri)));')
par(mar=c(1,1,3,1))
par(xpd=T)
plot.new()
plot(topo1,main="topo1")
plot(topo2,main="topo2")
plot(topo3,main="topo3")
par(mar=c(5,5,1,1))
boxplot(twisst4[,c(7:9)],xlab="Topology",ylab="Average weighting")
dev.off()


pdf("Twisst_MiashalleriKowaar_weightings_w100_paper.pdf",width=14,height=12,paper="special",pointsize=28)

par(mgp=c(3,0.75,0))
par(oma=c(2,2,2,2))
layout(matrix(c(1,2,3,4,5,5,5,5,5,5,5,5), nrow=3, byrow=TRUE),widths=c(0.25,1,1,1,1))
#topo1 (A_ly,(Mias_arenosa,(halleri_halleri,Kowa_arenosa)));
#topo2 (A_ly,(halleri_halleri,(Mias_arenosa,Kowa_arenosa)));
#topo3 (A_ly,(Kowa_arenosa,(Mias_arenosa,halleri)));

topo1 <-ape::read.tree(text='(A l,(Aa M,(Ah,Aa NM)));')
topo2<- ape::read.tree(text='(A l,(Ah,(Aa M,Aa NM)));')
topo3<- ape::read.tree(text='(A l,(Aa NM,(Aa M,Ah)));')
par(mar=c(1,1,3,1))
par(xpd=T)
plot.new()
plot(topo1,main="",font=c(3,2,3,2),cex=1.5,cex.lab=1.5,cex.main=1.5,font.main=1)
plot(topo2,main="",font=c(3,3,2,2),cex=1.5,cex.lab=1.5,cex.main=1.5,font.main=1)
plot(topo3,main="",font=c(3,2,2,3),cex=1.5,cex.lab=1.5,cex.main=1.5,font.main=1)
par(mar=c(5,5,1,1))
boxplot(twisst4[,c(7:9)],xlim=c(0.6,3.4),xlab="",ylab="Weighting",col="white",border=c("blue","black","red"),names=c("Topology 1","Topology 2","Topology 3"),cex=1.5,cex.axis=1.5,cex.lab=1.8,las=1)

dev.off()

pdf("Twisst_MiashalleriKowaar_weightings_w100_paper_top.pdf",width=14,height=4,paper="special",pointsize=28)

topo1 <-ape::read.tree(text='(A l,(Aa M,(Ah,Aa NM)));')
topo2<- ape::read.tree(text='(A l,(Ah,(Aa M,Aa NM)));')
topo3<- ape::read.tree(text='(A l,(Aa NM,(Aa M,Ah)));')
par(mar=c(1,1,3,1))
par(xpd=T)
plot.new()
par(mgp=c(3,0.75,0))
par(oma=c(2,2,2,2))
par(mfrow=c(1,3))
plot(topo1,main="",font=c(3,2,3,2),cex=1.5,cex.lab=1.5,cex.main=1.5,font.main=1)
plot(topo2,main="",font=c(3,3,2,2),cex=1.5,cex.lab=1.5,cex.main=1.5,font.main=1)
plot(topo3,main="",font=c(3,2,2,3),cex=1.5,cex.lab=1.5,cex.main=1.5,font.main=1)
dev.off()




require(qqman)
trace(manhattan,edit=T)
#change xaxis ticks:ticks <- tapply(d$pos, d$index, function(x) {min(x) + ((max(x) - min(x))/2)})

pdf("Twisst_MiashalleriKowaar_weightings_w100_Manhattanplot.pdf",width=10,height=10,paper="special",pointsize=20)
par(mar=c(5,5,1,1))
par(lwd=3)
par(mgp=c(2,0.75,0))
twisst_manhattanpl<-twisst4
twisst_manhattanpl$scaffold<-as.numeric(gsub("scaffold_","",twisst_manhattanpl$scaffold))
twisst_manhattanpl$snp<-paste0(twisst_manhattanpl$scaffold,twisst_manhattanpl$mid)
require(qqman)
manhattan(twisst_manhattanpl,snp="snp",chr="scaffold",bp="mid",p="topo3",ylab="Weighting",col=c("blue4","orange3"),suggestiveline=F,genomewideline=0.7,xlab="Scaffold",logp=F,ylim=c(0,1.1),cex.lab=1.5)

dev.off()


pdf("Twisst_MiashalleriKowaar_weightings_w100_Manhattanplot_per_scaffold.pdf",width=15,height=10,paper="special",pointsize=20)
par(mar=c(5,5,1,1))
par(lwd=3)
par(mgp=c(2,0.75,0))
for (i in 1:8)
	{twisst_manhattanpl<-twisst4[twisst4$scaffold==paste0("scaffold_",i),]
	twisst_manhattanpl$scaffold<-as.numeric(gsub("scaffold_","",twisst_manhattanpl$scaffold))
	twisst_manhattanpl$snp<-paste0(twisst_manhattanpl$scaffold,twisst_manhattanpl$mid)
	manhattan(twisst_manhattanpl,snp="snp",chr="scaffold",bp="mid",p="topo3",ylab="Weighting",col=c("blue4","orange3"),suggestiveline=F,genomewideline=0.7,xlab="Scaffold",logp=F,ylim=c(0,1.1),cex.lab=1.5)
	}
dev.off()

#HMA4
#scaffold_3	23475941	23483785
#MTP1
#scaffold_4	22776960	22778882
#HMA3
#scaffold_7	4960535	4963880
#NRAMP3
#AL4G13010	scaffold_4	1634727	1637073
#ZIP6<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=13949254&All3$Pos<=13950804,]
#ZIP11<-All3[All3$Scaffold=="scaffold_1"&All3$Pos>=29573151&All3$Pos<=29574444,]

End_pos_scaff<-tapply(twisst_manhattanpl$mid, twisst_manhattanpl$scaffold, function(x) {max(x)})
Scaffolds<-c(3,4,7,4,4,1)
Starts<-c(23475941,22776960,4960535,1634727,13949254,29573151)
Ends<-c(23483785,22778882,4963880,1637073,13950804,29574444)

pdf("Twisst_MiashalleriKowaar_weightings_w100_Manhattanplot_highlights.pdf",width=10,height=10,paper="special",pointsize=20)
par(mar=c(5,5,1,1))
par(lwd=3)
par(mgp=c(2,0.75,0))
twisst_manhattanpl<-twisst4
twisst_manhattanpl$scaffold<-as.numeric(gsub("scaffold_","",twisst_manhattanpl$scaffold))
twisst_manhattanpl$snp<-paste0(twisst_manhattanpl$scaffold,twisst_manhattanpl$mid)
require(qqman)
manhattan(twisst_manhattanpl,snp="snp",chr="scaffold",bp="mid",p="topo3",ylab="Weighting",col=c("blue4","orange3"),suggestiveline=F,genomewideline=0.7,xlab="Scaffold",logp=F,ylim=c(0,1.1),cex.lab=1.5)
for(j in 1:5)
	{scaffold=Scaffolds[j]
	start=Starts[j]
	end=Ends[j]
	points((sum(End_pos_scaff[1:(scaffold-1)])+twisst_manhattanpl$mid[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end]),
	twisst_manhattanpl$topo3[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end],col="red",pch=19)
	}
scaffold=Scaffolds[6]
start=Starts[6]
end=Ends[6]
points((twisst_manhattanpl$mid[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end]),
twisst_manhattanpl$topo3[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end],col="red",pch=19)

dev.off()

pdf("Twisst_MiashalleriKowaar_weightings_w100_Manhattanplot_highlights_per_scaffold.pdf",width=15,height=10,paper="special",pointsize=20)
par(mar=c(5,5,1,1))
par(lwd=3)
par(mgp=c(3,0.75,0))
for (i in 1:8)
	{twisst_manhattanpl<-twisst4[twisst4$scaffold==paste0("scaffold_",i),]
	twisst_manhattanpl$scaffold<-as.numeric(gsub("scaffold_","",twisst_manhattanpl$scaffold))
	twisst_manhattanpl$snp<-paste0(twisst_manhattanpl$scaffold,twisst_manhattanpl$mid)
	manhattan(twisst_manhattanpl,snp="snp",chr="scaffold",bp="mid",p="topo3",ylab="Weighting",col=c("blue4","orange3"),suggestiveline=F,genomewideline=0.7,xlab="Scaffold (Mb)",logp=F,ylim=c(0,1.1),cex.lab=1.5,xaxt="n",cex.axis=1.3)
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="g"),cex=1.5,cex.lab=1.5,cex.axis=1.3)
	for(j in 1:6)
		{scaffold=Scaffolds[j]
		start=Starts[j]
		end=Ends[j]
		if (scaffold==i)
			{points(twisst_manhattanpl$mid[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end],
			twisst_manhattanpl$topo3[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end],col="red",pch=19)
			}
		}
	box()
	}

dev.off()

pdf("Twisst_MiashalleriKowaar_weightings_w100_Manhattanplot_highlights_per_scaffold_lines.pdf",width=15,height=10,paper="special",pointsize=20)
par(mar=c(5,5,3,1))
par(lwd=3)
par(mgp=c(3,0.75,0))
for (i in 1:8)
	{twisst_manhattanpl<-twisst4[twisst4$scaffold==paste0("scaffold_",i),]
	twisst_manhattanpl$scaffold<-as.numeric(gsub("scaffold_","",twisst_manhattanpl$scaffold))
	twisst_manhattanpl$snp<-paste0(twisst_manhattanpl$scaffold,twisst_manhattanpl$mid)
	manhattan(twisst_manhattanpl,snp="snp",chr="scaffold",bp="mid",p="topo3",ylab="Weighting",col=c("blue4","orange3"),suggestiveline=F,genomewideline=0.7,xlab="Scaffold (Mb)",logp=F,ylim=c(0,1.1),cex.lab=1.5,xaxt="n",cex.axis=1.3,type="l",lwd=3,main=paste0("Scaffold ",i))
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="g"),cex=1.5,cex.lab=1.5,cex.axis=1.3)
	for(j in 1:6)
		{scaffold=Scaffolds[j]
		start=Starts[j]
		end=Ends[j]
		if (scaffold==i)
			{points(twisst_manhattanpl$mid[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end],
			twisst_manhattanpl$topo3[twisst_manhattanpl$scaffold==scaffold&twisst_manhattanpl$end>start&twisst_manhattanpl$start<end],col="red",pch=19,cex=0.7)
			}
		}
	box()
	}

dev.off()

#14470836	14471310
par(mar=c(5,5,3,1))
par(lwd=3)
par(mgp=c(3,0.75,0))
for (i in 3)
	{twisst_manhattanpl<-twisst4[twisst4$scaffold==paste0("scaffold_",i),]
	twisst_manhattanpl$scaffold<-as.numeric(gsub("scaffold_","",twisst_manhattanpl$scaffold))
	twisst_manhattanpl$snp<-paste0(twisst_manhattanpl$scaffold,twisst_manhattanpl$mid)
	manhattan(twisst_manhattanpl,snp="snp",chr="scaffold",bp="mid",p="topo3",ylab="Weighting",col=c("blue4","orange3"),suggestiveline=F,genomewideline=0.7,xlab="Scaffold (Mb)",logp=F,ylim=c(0,1.1),cex.lab=1.5,xaxt="n",cex.axis=1.3,type="l",lwd=3,main=paste0("Scaffold ",i))
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="g"),cex=1.5,cex.lab=1.5,cex.axis=1.3)
	points(twisst_manhattanpl$mid[twisst_manhattanpl$scaffold=="scaffold_3"&twisst_manhattanpl$end>14470836&twisst_manhattanpl$start<14473310],
	twisst_manhattanpl$topo3[twisst_manhattanpl$scaffold=="scaffold_3"&twisst_manhattanpl$end>14470836&twisst_manhattanpl$start<14473310],col="red",pch=19,cex=0.7)
	box()
	}

#sign.
nrow(twisst4[twisst4$topo3>=0.6,])
#1273
nrow(twisst4[twisst4$topo1>=0.6,])
#1942

nrow(twisst4[twisst4$topo3>=0.7,])
#433
nrow(twisst4[twisst4$topo1>=0.7,])
#537

Marhalleri_05<-twisst4[twisst4$topo3>=0.7,]
nrow(Marhalleri_05)/nrow(twisst4)*100
#0.716662

Karhalleri_05<-twisst4[twisst4$topo1>=0.7,]

pdf("Hist_Twisst_MiashalleriKowaar_topo3_w100_paper.pdf",width=15,height=8,paper="special",pointsize=20)
par(oma=c(1,1,1,1))
par(mgp=c(2.5,0.75,0))
plot(hist(twisst4$topo3),main="Introgression",cex=1.5,cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="Topology 3 weighting",las=1,cex.axis=1.2,yaxt="n",ylab="",col="black")
mgp.axis(side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = "."),cex=1,cex.lab=1.2,cex.axis=1.2,las=1)
mtext(side=2,"Frequency",line=-0.5,outer=T,cex=1.5)
abline(v=0.7,col="red")
dev.off()

require(Hmisc)
pdf("Hist_Twisst_MiashalleriKowaar_topo3_w100_paper_square.pdf",width=10,height=8,paper="special",pointsize=30)
par(oma=c(1,1,1,1))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
hist(twisst4$topo3,breaks=100,main="",cex.lab=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,4000))
mgp.axis(pos=min(hist(twisst4$topo3,breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2)/1000,big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.2,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=axTicks(side=1),cex=1,cex.lab=1.2,cex.axis=1.2,las=1)
mtext(side=2,"Frequency",line=-0.5,outer=T,cex=1.5)
mtext(side=1,"Twisst weighting (Topology 3)",line=-2.5,outer=T,cex=1.5)
lines(c(0.7,0.7),c(0,4000),col="red",xpd=F,lwd=3)
dev.off()



###################################################
#find overlaps between windows and A. lyrata genes#
###################################################
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

### Step 1: import genes from GFF files

lyr_txdm = makeTxDbFromGFF(file ="D:/Lyrata/Alyrata_384_v2.1.gene.gff3", format = c("gff3"), organism = "Arabidopsis lyrata")

#lyr_txdm
lyrGenes=genes(lyr_txdm)
names(lyrGenes)<-NULL

#Define a function that allows the addition of promotor regions
extend <- function(x, upstream=0, downstream=0) #from https://support.bioconductor.org/p/78652/
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

# Use the function to extend promotor region 2kb
  # Upstream=2000 adds 2kb upstream from transcription start site
  # Downstream=0 adds 0bp downstream from the gene.
LyrGenes_plusPromot = extend(lyrGenes,upstream=2000, downstream=0)

### Convert list of significant markers to a "Genomic Ranges" table. 
Twisst_sig_Marhalleri_05_GRange<-GRanges(seqnames=tolower(Marhalleri_05$scaffold),ranges=IRanges(start=Marhalleri_05$start,end=Marhalleri_05$end))
values(Twisst_sig_Marhalleri_05_GRange)<-Marhalleri_05[,4:9]
Twisst_sig_Karhalleri_05_GRange<-GRanges(seqnames=tolower(Karhalleri_05$scaffold),ranges=IRanges(start=Karhalleri_05$start,end=Karhalleri_05$end))
values(Twisst_sig_Karhalleri_05_GRange)<-Karhalleri_05[,4:9]

Twisst_Marhalleri<-findOverlaps(Twisst_sig_Marhalleri_05_GRange,LyrGenes_plusPromot)
Twisst_Marhalleri2<-pintersect(Twisst_sig_Marhalleri_05_GRange[queryHits(Twisst_Marhalleri)],LyrGenes_plusPromot[subjectHits(Twisst_Marhalleri)])
Twisst_Marhalleri3<-as.data.frame(Twisst_Marhalleri2)
Twisst_Marhalleri_lyrGenes=mergeByOverlaps(Twisst_sig_Marhalleri_05_GRange,LyrGenes_plusPromot,type=c("any"))
Twisst_Marhalleri_lyrGenes_df=data.frame(as.data.frame(Twisst_Marhalleri_lyrGenes$Twisst_sig_Marhalleri_05_GRange),as.data.frame(Twisst_Marhalleri_lyrGenes$LyrGenes_plusPromot),Twisst_Marhalleri3$width)
colnames(Twisst_Marhalleri_lyrGenes_df)=c("Scaffold", "Start", "End", "Width", "Strand", "Middle","Num_SNPs","lnL","topo1","topo2","topo3","Scaffold_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap")

Twisst_Karhalleri<-findOverlaps(Twisst_sig_Karhalleri_05_GRange,LyrGenes_plusPromot)
Twisst_Karhalleri2<-pintersect(Twisst_sig_Karhalleri_05_GRange[queryHits(Twisst_Karhalleri)],LyrGenes_plusPromot[subjectHits(Twisst_Karhalleri)])
Twisst_Karhalleri3<-as.data.frame(Twisst_Karhalleri2)
Twisst_Karhalleri_lyrGenes=mergeByOverlaps(Twisst_sig_Karhalleri_05_GRange,LyrGenes_plusPromot,type=c("any"))
Twisst_Karhalleri_lyrGenes_df=data.frame(as.data.frame(Twisst_Karhalleri_lyrGenes$Twisst_sig_Karhalleri_05_GRange),as.data.frame(Twisst_Karhalleri_lyrGenes$LyrGenes_plusPromot),Twisst_Karhalleri3$width)
colnames(Twisst_Karhalleri_lyrGenes_df)=c("Scaffold", "Start", "End", "Width", "Strand", "Middle","Num_SNPs","lnL","topo1","topo2","topo3","Scaffold_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap")

require(openxlsx)
Thal_description<-read.xlsx("E:/Lyrata_RBH/Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx")
UtesMetal=read.xlsx("D:/Lyrata/metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6:16)]

Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(Twisst_Marhalleri_lyrGenes_df,Thal_description, by.x="Gene", by.y="Alyr_ID",all.x=TRUE)
colnames(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup=merge(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup<-Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup[,c(2:13,1,14:40)]

write.xlsx(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_Mias_Twisst_w100_MahalleriKaar_07.xlsx",overwrite=T)


Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup_Thalgenes= merge(Twisst_Karhalleri_lyrGenes_df,Thal_description, by.x="Gene", by.y="Alyr_ID",all.x=TRUE)
colnames(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup_Thalgenes)[1]="Lyr_Gene"

Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup=merge(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup_Thalgenes,UtesMetal, by.x ="Ath_ID", by.y="AGI_Number",all.x=TRUE)
Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup<-Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup[,c(2:13,1,14:40)]

write.xlsx(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup,"GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_Kowa_Twisst_w100_MahalleriKaar_07.xlsx",overwrite=T)



length(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup$Lyr_Gene[!duplicated(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup$Lyr_Gene)])
#565
length(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup$Lyr_Gene[!duplicated(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup$Lyr_Gene)])
#811

length(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup$Name[!duplicated(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup$Name)])
#11
length(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup$Name[!duplicated(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup$Name)])
#17

paste0(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup$Name[!duplicated(Twisst_Marhalleri_lyrGenes_df_desc_w_orthogroup$Name)],sep=",",collapse="")
# "NA,ZIP11,SERAT2;1,PIC1,HMA4,NRAMP3,ZIP6,MTP1,SERAT2;2,HMA2,HMA3,"

paste0(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup$Name[!duplicated(Twisst_Karhalleri_lyrGenes_df_desc_w_orthogroup$Name)],sep=",",collapse="")
#"NA,CAX11,FRD3,GufA,ZTP50,PCR5,PCR4,PCR7,NRAMP5,IRT2,CYP82C4,ZIFl1,CAX8,CAX7,FRO5,RAN1,bHLH105, ILR3,"







genelist<-read.table("F:/Cu_project/Lyrata/Alyrata_384_v2.1.gene.gff3",sep="\t",header=F)
colnames(genelist)<-c("Scaffold","Source","Type","gene_start","gene_end","Score","Strand","Phase","Attributes")
genelist<-genelist[genelist$Type=="gene",]
ID<-substr(as.character(genelist$Attributes),4,nchar(as.character(genelist$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
genelist$Lyr_Gene<-ID3[,1]

genelist2<-read.table("F:/Cu_project/Lyrata/Alyrata_384_v2.1.gene_exons.gff3",sep="\t",header=F)
colnames(genelist2)<-c("Scaffold","Source","Type","gene_start","gene_end","Score","Strand","Phase","Attributes")
exonlist<-na.exclude(genelist2[genelist2$Type=="exon",])
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
ID<-substr(as.character(exonlist$Attributes),4,nchar(as.character(exonlist$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
exonlist$Lyr_Gene<-ID3[,1]
exonlist<-exonlist[!duplicated(exonlist[,-9]),]


pdf("Twisst_Mias_weightings_w100_relative_HMA4_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_3"
Start<-23475941
End<-23483785
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G52820",]

par(lwd=3)
par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-1#-

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()


#HMA2
#scaffold_7	4971705	4976723
pdf("Twisst_Mias_weightings_w100_relative_HMA2_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_7"
Start<-4971705
End<-4976723
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G22090",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-2#+

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()


#HMA3
#scaffold_7	4960535	4963880
pdf("Twisst_Mias_weightings_w100_relative_HMA3_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_7"
Start<-4960535
End<-4963880
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G22080",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-2#+

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#MTP1
#scaffold_4	22776960	22778882
pdf("Twisst_Mias_weightings_w100_relative_MTP1_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_4"
Start<-22776960
End<-22778882
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL4G46270",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-2#+

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()



#ZIP6 AL4G24850
#scaffold_4	phytozomev11	gene	13949254	13950804
pdf("Twisst_Mias_weightings_w100_relative_ZIP6_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_4"
Start<-13949254
End<-13950804
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL4G24850",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-1#-

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#ZIP11 AL1G64050
#AT1G55910	scaffold_1	phytozomev11	gene	29573151	29574444
pdf("Twisst_Mias_weightings_w100_relative_ZIP11_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_1"
Start<-29573151
End<-29574444
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL1G64050",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-2#+

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#NRAMP3 AL4G13010
#scaffold_4	phytozomev11	gene	1634727	1637073
pdf("Twisst_Mias_weightings_w100_relative_NRAMP3_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_4"
Start<-1634727
End<-1637073
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL4G13010",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-1#-

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#IRT2 AL7G34360 scaffold_7	phytozomev11	gene	10093506	10094856
pdf("Twisst_Mias_weightings_w100_relative_IRT2_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_7"
Start<-10093506
End<-10094856
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G34360",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-1#-

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-100000))&(twisst4$mid<(End+windowsize+100000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#ZIP8 AL8G14120 scaffold_8	phytozomev11	gene	2517313	2519529

pdf("Twisst_Mias_weightings_w100_relative_ZIP8_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_8"
Start<-2517313
End<-2519529
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL8G14120",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-2#+

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#FRO5 AL6G35310 scaffold_6	phytozomev11	gene	10296305	10299581

pdf("Twisst_Mias_weightings_w100_relative_FRO5_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_6"
Start<-10296305
End<-10299581
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL6G35310",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-1#-

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#FER1 AL6G10450 scaffold_6	phytozomev11	gene	187631	189383

pdf("Twisst_Mias_weightings_w100_relative_FER1_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_6"
Start<-187631
End<-189383
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL6G10450",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-2#+

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()

#FRD3 AL3G19430 scaffold_3	phytozomev11	gene	3352004	3357699

pdf("Twisst_Mias_weightings_w100_relative_FRD3_paper2.pdf",width=12,height=7,paper="special",pointsize=19)
Scaff<-"scaffold_3"
Start<-3352004
End<-3357699
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G19430",]

par(mar=c(6,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,1,0))
arrowdir2<-1#-

twisst4_scaff7<-twisst4[(twisst4$scaffold==Scaff)&(twisst4$mid>(Start-windowsize-10000))&(twisst4$mid<(End+windowsize+10000)),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="",ylab="Weighting",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,yaxt="n",xlim=c(Start-windowsize,End+windowsize))
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3)
mgp.axis(side=2,at=axTicks(side=2)[-length(axTicks(side=2))],labels=formatC(axTicks(side=2)[-length(axTicks(side=2))],format="f",decimal.mark = ".",digits=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue",lwd=3)
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black",lwd=3)
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red",lwd=3)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=3)
	}
arrows(x0=Start,y0=1.1,x1=End,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=3)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((End-Start)/2),1.2,"HMA4",col="black",font=3)
#legend("topleft",lwd=2,col=c("blue","black","red"),legend=c("Topology 1","Topology 2","Topology 3"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)
box(,lwd=1)
dev.off()



source("plot_twisst.R")
cols = c(
"#0075DC", #Blue
"#2BCE48", #Green
"#FFA405", #Orpiment
"#5EF1F2", #Sky
"#FF5005", #Zinnia
"#005C31", #Forest
"#00998F", #Turquoise
"#FF0010", #Red
"#9DCC00", #Lime
"#003380", #Navy
"#F0A3FF", #Amethyst
"#740AFF", #Violet
"#426600", #Quagmire
"#C20088", #Mallow
"#94FFB5") #Jade
par(mfrow =c(1,1), mar = c(3,3,1,1), xpd = FALSE)
plot.weights(weights_dataframe=weights_by_chrom4[[3]],positions=cbind(window_data_by_chrom4[[3]]$start,window_data_by_chrom4[[3]]$end),
  line_cols=NA, lwd= 0, fill_cols=cols, xlim =c((23475941-30000),(23485785+30000)),stacked=TRUE,ylim=c(0,1.2))

pdf("Twisst_Mias_weightings_w25_relative_HMA4_Smoothed_plot.pdf",width=30,height=20,paper="special",pointsize=20)
weights_smooth = smooth.weights(window_positions=window_data_by_chrom4[[3]]$mid,weights_dataframe = weights_by_chrom4[[3]],
  span = 0.001, window_sites=window_data_by_chrom4[[3]]$sites)
plot.weights(weights_dataframe=weights_smooth,positions=cbind(window_data_by_chrom4[[3]]$start,window_data_by_chrom4[[3]]$end),
  line_cols=NA, lwd= 0, fill_cols=cols, xlim =c((23475941-30000),(23485785+30000)),stacked=TRUE,ylim=c(0,1.2))
arrowdir2<-1#-
arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1.1,x1=Lyr_Genesofinterest$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=23475941,y0=1.1,x1=23485785,y1=1.1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1.1,x1=Exonsofinterest$gene_end[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1.1,x1=Exonsofinterestcand$gene_end[i],y1=1.1,code=arrowdir2,length=0,col="red",lwd=3)
	}
dev.off()




#MT1c
#AL1G17780	scaffold_1	2785797	2788520
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_MT1c.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_1"
Start<-2785797
End<-2788520
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL1G17780",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()

#NRAMP3
#AL4G13010	scaffold_4	1634727	1639073
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_NRAMP3.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_4"
Start<-1634727
End<-1639073
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL4G13010",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()

#bHLH011
#AL7G14990	scaffold_7	2048538	2052007
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_bHLH011.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_7"
Start<-2048538
End<-2052007
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G14990",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()

#FER1
#AL6G10450	scaffold_6	185631	189383
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_FER1.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_6"
Start<-185631
End<-189383
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL6G10450",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()

#ZIP8
#AL8G14120	scaffold_8	2515313	2519529
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_ZIP8.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_8"
Start<-2515313
End<-2519529
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL8G14120",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()


#HTR4
#AL6G41640	scaffold_6	17880366	17883252
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_HTR4.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_6"
Start<-17880366
End<-17883252
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL6G41640",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()


#FRO2
#AL1G10810 scaffold_1	266658	272551
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_FRO2.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_1"
Start<-266658
End<-272551
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL1G10810",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()

#PIC1
#AL3G47010 scaffold_3	19789770	19793689
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_PIC1.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_3"
Start<-19789770
End<-19793689
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G47010",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()

#MEB3
#AL7G24770 scaffold_7	5991780	5998124
pdf("Twisst_MiashalleriKowaarth_weightings_w100_relative_MEB3.pdf",width=30,height=20,paper="special",pointsize=20)
Scaff<-"scaffold_7"
Start<-5991780
End<-5998124
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
Lyr_Genes1<-genelist[genelist$Scaffold==Scaff,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(Start-windowsize)&Lyr_Genes1$gene_start<=(Start+windowsize)|Lyr_Genes1$gene_end>=(End-windowsize)&Lyr_Genes1$gene_end<=(End+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==Scaff,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(Start-windowsize)&Exonsofinterest1$gene_start<=(Start+windowsize)|Exonsofinterest1$gene_end>=(End-windowsize)&Exonsofinterest1$gene_end<=(End+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G24770",]

arrowdir2<-1#-

twisst4_scaff7<-twisst4[twisst4$scaffold==Scaff&twisst4$mid<(Start+windowsize)&twisst4$mid>(End-windowsize),]
plot(twisst4_scaff7$topo1~twisst4_scaff7$mid,col="blue",type="l",ylim=c(0,1.2),xlab="Scaffold position",ylab="Absolute weight")
rect(Start,-1,End,1.5,col="lightgrey",border = NA)	
lines(twisst4_scaff7[,7]~twisst4_scaff7$mid,col="blue")
lines(twisst4_scaff7[,8]~twisst4_scaff7$mid,col="black")
lines(twisst4_scaff7[,9]~twisst4_scaff7$mid,col="red")

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=1,x1=Lyr_Genesofinterest$gene_end[i],y1=1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
	}
arrows(x0=Start,y0=1,x1=End,y1=1,code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=1,x1=Exonsofinterest$gene_end[i],y1=1,code=1,length=0,col="grey",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=1,x1=Exonsofinterestcand$gene_end[i],y1=1,code=arrowdir2,length=0,col="red",lwd=5)
	}
dev.off()












dist_list<-paste0("Mias.w25.",chromNames,".dist")
dist_by_chrom4 = list()
for (i in 1:length(dist_list)){
    dist <- read.table(dist_list[i], header = T, as.is = T)
    dist_by_chrom4[[i]] <- dist
     }
dist_all4=data.frame()
for(i in 1:8){
dist_all4=rbind(dist_all4,dist_by_chrom4[[i]])}
dist_all4<-dist_all4[rownames(weights_all4),]

summary(dist_all4[,grep("Mias_arenosa_Mias_halleri",colnames(dist_all4))])

plot(dist_all4[,grep("Mias_arenosa_Mias_halleri",colnames(dist_all4))[1]]~twisst4$mid,col="black",type="l",ylim=c(0,2),xlab="Scaffold position",ylab="Distance",xlim=c((23475941-30000),(23485785+30000)))
rect(23475941,0,23485785,2,col="lightgrey",border = NA)	
for (i in (grep("Mias_arenosa_Mias_halleri",colnames(dist_all4),invert=T)))
	{lines(dist_all4[,i]~twisst4$mid,col="black")
	}
for (i in (grep("Mias_arenosa_Mias_halleri",colnames(dist_all4))))
	{lines(dist_all4[,i]~twisst4$mid,col="red")
	}










twisst4['data']='4c'
twisst4['source']='All'
twisst4$start=as.numeric(twisst4$start)
twisst4$end=as.numeric(twisst4$end)
```


```{r}
#Define locations of meiosis outliers
topoNames = names(weights)
regions <- matrix(c("scaffold_1",0,0.5e6,220193,225616,
                    "scaffold_1",9.5e6,10.5e6,9779387,9791542,
                    "scaffold_2",17e6,18e6,17706052,17715900,
                    "scaffold_2",12e6,13e6,12427364,12431683,
                    "scaffold_4",11e6,12e6,11123922,11131646,
                    "scaffold_4",22.5e6,23.5e6,22845898,22850613,
                    "scaffold_6", 1.75e6,2.25e6,2001440,2005979
                    ), ncol=5, byrow=T)
reg_names=c("PRD3","ZYP1a/b","PDS5","ASY1","SMC3","ASY3","SYN1")
plot_order = c(1,2,3)
# Get the weights for the windows within the meoisis genes
m_winds3=data.frame()
for (x in 1:nrow(regions)){
    df = twisst3[twisst3$scaffold==regions[x,1] & twisst3$start < as.numeric(regions[x,4]) & twisst3$end > as.numeric(regions[x,5]),]
    m_winds3 = rbind(m_winds3,df)
}
m_winds3=rbind(m_winds3,twisst3[twisst3$scaffold=="scaffold_2" & twisst3$topo1 == weights_by_chrom3[[2]]["447",1],]) 
m_winds4=data.frame()
for (x in 1:nrow(regions)){
    df = twisst4[twisst4$scaffold==regions[x,1] & twisst4$start < as.numeric(regions[x,4]) & twisst4$end > as.numeric(regions[x,5]),]
    m_winds4 = rbind(m_winds4,df)
}
m_winds4=rbind(m_winds4,twisst4[twisst4$scaffold=="scaffold_2" & twisst4$topo1 == weights_by_chrom4[[2]]["447",1],]) 
m_winds4['data']='4c'
m_winds3['data']='3c'
m_winds4['source']='Mei'
m_winds3['source']='Mei'
#Make master data frame with meosis weights along with all other window weights
twisst3=rbind(twisst3,m_winds3)
twisst4=rbind(twisst4,m_winds4)
twisst=rbind(twisst3,twisst4)
twisst[twisst$data=="4c",]$data="Baltic"
twisst[twisst$data=="3c",]$data="S. Carpathians"
m.twisst=melt(twisst[,c("source","data","topo1","topo2","topo3")],id.vars=c("source","data"))
#Plotting command
ggplot(m.twisst,aes(x=variable,y=value,fill=source))+geom_boxplot(outlier.size=0.5)+theme_bw()+ylab("Weights")+facet_grid(~data)+theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1,size=14),legend.title=element_blank(),legend.text=element_text(size=14),axis.title.y=element_text(size=14),strip.text.x=element_text(size=14),axis.text.y=element_text(size=12))
```

# Weights are normally distributed for the most part.  
```{r}
ggplot(m.twisst, aes(y = data, fill = variable, x = value)) + geom_density_ridges(alpha = 0.5) + facet_grid(~source)
```

# Significant differences between groups in above plot for both parametric and nonparametric tests
```{r}
m.twisst %>% group_by(data, variable) %>% summarize(p.t = t.test(value ~ source)$p.value, stat.t = t.test(value ~ source)$statistic, p.w = wilcox.test(value ~ source)$p.value, stat.w = wilcox.test(value ~ source)$statistic)
```














for (i in 1:8)
	{assign(paste0("trees_scaff_",i),read.tree(paste0("MiashalleriKowaarenosa.topologies.w100.m1.scaffold_",i,".trees")))
	}
ggdensitree(trees_scaff_1)


MhaKar_trees_scaff1<-read.tree("MiashalleriKowaarenosalyr.Vcf.scfscaffold_1.w100m1.LDphase.phyml_bionj.trees")





for (i in 1:length(twisst_data$topo1)) twisst_data$topo1[[i]] <- ladderize(unroot(twisst4$topo1[[i]]))

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "unrooted", edge.color=topo_cols[n], edge.width=5, rotate.tree = 90, cex = 1, adj = .5, label.offset=.2)
  mtext(side=3,text=paste0("topo",n))
  }

#rooted topologies
source("../plot_twisst_windows_LS.R")
twisst_data <- import.twisst(weights_files=weights_files4,window_data_files=window_data_files,names=chromNames)
for (i in 1:length(twisst_data$topos)) twisst_data$topos[[i]] <- root(twisst_data$topos[[i]], "A_ly", resolve.root = T)

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "clad", edge.color=topo_cols[n], edge.width=5, label.offset=.1, cex = 1)
  mtext(side=3,text=paste0("topo",n))
  }




tree1<-read.tree("MiashalleriKowaarenosalyr.Vcf.scfscaffold_1.w100m1.LDphase.phyml_bionj.no_NAs.txt")
names(tree1)<-paste0("window_",i:length(tree1))


tree1<-phytools::read.newick("MiashalleriKowaarenosalyr.Vcf.scfscaffold_1.w100m1.LDphase.phyml_bionj.trees.gz")

require(ggtree)
#from precrec
df <- fortify(tree1)
p+geom_tree(data=df, color='firebrick')


ggdensitree(tree1)
trace(?, edit=TRUE)
p <- ggtree(tree1,layout="rectangular",color="lightblue",alpha=.3)
pdf("Twisst_MiashalleriKowaarenosa_w100_tree_scaff1.pdf",width=50,height=30,paper="special",pointsize=20)
plot(p)
dev.off()


ape::write.tree(p, file='Twisst_MiashalleriKowaarenosa_w100_tree_scaff1.txt')
ggsave('Twisst_MiashalleriKowaarenosa_w100_tree_scaff1_ggtree.pdf', width = 50, height = 50, units = "cm", limitsize = FALSE)

The resulting 18,000 trees were combined to a maximum clade credibility tree in TreeAnnotator66 version 1.7.5 and visualized in FigTree66 version 1.4.1.







