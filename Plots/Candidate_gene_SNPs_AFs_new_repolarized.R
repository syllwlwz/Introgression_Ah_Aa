Buko_arenosa<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Buko_arenosa_Treemix_all.vcf",sep="\t",header=T)
Kato_arenosa<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Kato_arenosa_Treemix_all.vcf",sep="\t",header=T)
Mias_arenosa<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Mias_arenosa_Treemix_all.vcf",sep="\t",header=T)
Piek_arenosa<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Piek_arenosa_Treemix_all.vcf",sep="\t",header=T)
Kowa_arenosa<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Kowa_arenosa_Treemix_all.vcf",sep="\t",header=T)
Zapa_arenosa<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Zapa_arenosa_Treemix_all.vcf",sep="\t",header=T)

Buko_halleri<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Buko_halleri_Treemix_all.vcf",sep="\t",header=T)
Mias_halleri<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Mias_halleri_Treemix_all.vcf",sep="\t",header=T)
Piek_halleri<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Piek_halleri_Treemix_all.vcf",sep="\t",header=T)
Kowa_halleri<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Kowa_halleri_Treemix_all.vcf",sep="\t",header=T)
Zapa_halleri<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/Zapa_halleri_Treemix_all.vcf",sep="\t",header=T)


thaliana<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/thaliana_Treemix_all.vcf",sep="\t",header=T)
lyrata<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/lyrata_Treemix_all.vcf",sep="\t",header=T)
cebennensis<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/cebennensis_all_Treemix.vcf",sep="\t",header=T)
croatica<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/croatica_all_Treemix.vcf",sep="\t",header=T)
pedemontana<-read.delim("E:/Introgression/Demography_redo/All5/Repolarized/pedemontana_Treemix_all.vcf",sep="\t",header=T)

thaliana<-thaliana[complete.cases(thaliana[,1:4]),1:4]
thaliana$AC<-as.numeric(thaliana$AC)
lyrata<-lyrata[complete.cases(lyrata[,1:4]),1:4]
lyrata$AC<-as.numeric(lyrata$AC)
cebennensis<-cebennensis[complete.cases(cebennensis[,1:4]),1:4]
cebennensis$AC<-as.numeric(cebennensis$AC)
croatica<-croatica[complete.cases(croatica[,1:4]),1:4]
croatica$AC<-as.numeric(croatica$AC)
pedemontana<-pedemontana[complete.cases(pedemontana[,1:4]),1:4]
pedemontana$AC<-as.numeric(pedemontana$AC)

Buko_arenosa<-Buko_arenosa[complete.cases(Buko_arenosa[,1:4]),1:4]
Buko_arenosa$AC<-as.numeric(Buko_arenosa$AC)
Kato_arenosa<-Kato_arenosa[complete.cases(Kato_arenosa[,1:4]),1:4]
Kato_arenosa$AC<-as.numeric(Kato_arenosa$AC)
Mias_arenosa<-Mias_arenosa[complete.cases(Mias_arenosa[,1:4]),1:4]
Mias_arenosa$AC<-as.numeric(Mias_arenosa$AC)
Piek_arenosa<-Piek_arenosa[complete.cases(Piek_arenosa[,1:4]),1:4]
Piek_arenosa$AC<-as.numeric(Piek_arenosa$AC)
Kowa_arenosa<-Kowa_arenosa[complete.cases(Kowa_arenosa[,1:4]),1:4]
Kowa_arenosa$AC<-as.numeric(Kowa_arenosa$AC)
Zapa_arenosa<-Zapa_arenosa[complete.cases(Zapa_arenosa[,1:4]),1:4]
Zapa_arenosa$AC<-as.numeric(Zapa_arenosa$AC)

Buko_halleri<-Buko_halleri[complete.cases(Buko_halleri[,1:4]),1:4]
Buko_halleri$AC<-as.numeric(Buko_halleri$AC)
Mias_halleri<-Mias_halleri[complete.cases(Mias_halleri[,1:4]),1:4]
Mias_halleri$AC<-as.numeric(Mias_halleri$AC)
Piek_halleri<-Piek_halleri[complete.cases(Piek_halleri[,1:4]),1:4]
Piek_halleri$AC<-as.numeric(Piek_halleri$AC)
Kowa_halleri<-Kowa_halleri[complete.cases(Kowa_halleri[,1:4]),1:4]
Kowa_halleri$AC<-as.numeric(Kowa_halleri$AC)
Zapa_halleri<-Zapa_halleri[complete.cases(Zapa_halleri[,1:4]),1:4]
Zapa_halleri$AC<-as.numeric(Zapa_halleri$AC)

Buko_arenosa$AC[Buko_arenosa$AN==0]<-0
Kato_arenosa$AC[Kato_arenosa$AN==0]<-0
Mias_arenosa$AC[Mias_arenosa$AN==0]<-0
Piek_arenosa$AC[Piek_arenosa$AN==0]<-0
Kowa_arenosa$AC[Kowa_arenosa$AN==0]<-0
Zapa_arenosa$AC[Zapa_arenosa$AN==0]<-0

Buko_halleri$AC[Buko_halleri$AN==0]<-0
Mias_halleri$AC[Mias_halleri$AN==0]<-0
Piek_halleri$AC[Piek_halleri$AN==0]<-0
Kowa_halleri$AC[Kowa_halleri$AN==0]<-0
Zapa_halleri$AC[Zapa_halleri$AN==0]<-0

thaliana$AN<-as.numeric(thaliana$AN)
lyrata$AN<-as.numeric(lyrata$AN)
cebennensis$AN<-as.numeric(cebennensis$AN)
croatica$AN<-as.numeric(croatica$AN)
pedemontana$AN<-as.numeric(pedemontana$AN)

All<-cbind(thaliana[,1:2],thaliana[,3]/thaliana[,4],lyrata[,3]/lyrata[,4],cebennensis[,3]/cebennensis[,4],pedemontana[,3]/pedemontana[,4],croatica[,3]/croatica[,4],
Buko_arenosa[,3]/Buko_arenosa[,4],Kato_arenosa[,3]/Kato_arenosa[,4],Mias_arenosa[,3]/Mias_arenosa[,4],
Piek_arenosa[,3]/Piek_arenosa[,4],Kowa_arenosa[,3]/Kowa_arenosa[,4],Zapa_arenosa[,3]/Zapa_arenosa[,4],Buko_halleri[,3]/Buko_halleri[,4],Mias_halleri[,3]/Mias_halleri[,4],
Piek_halleri[,3]/Piek_halleri[,4],Kowa_halleri[,3]/Kowa_halleri[,4],Zapa_halleri[,3]/Zapa_halleri[,4])

colnames(All)<-c("Scaffold","Pos","thaliana_AF","lyrata_AF","cebennensis_AF","pedemontana_AF","croatica_AF","Buko_arenosa_AF","Kato_arenosa_AF","Mias_arenosa_AF",
"Piek_arenosa_AF","Kowa_arenosa_AF","Zapa_arenosa_AF","Buko_halleri_AF","Mias_halleri_AF","Piek_halleri_AF","Kowa_halleri_AF","Zapa_halleri_AF")

All2<-All[!(apply(All[,8:18],1,sum)==0),]
#3,860,675 of 5,278,999
All3<-All2[complete.cases(All2),]
#3,859,234
require(RColorBrewer)


All_wooutgroup<-cbind(Buko_arenosa[,1:2],
Buko_arenosa[,3]/Buko_arenosa[,4],Kato_arenosa[,3]/Kato_arenosa[,4],Mias_arenosa[,3]/Mias_arenosa[,4],
Piek_arenosa[,3]/Piek_arenosa[,4],Kowa_arenosa[,3]/Kowa_arenosa[,4],Zapa_arenosa[,3]/Zapa_arenosa[,4],Buko_halleri[,3]/Buko_halleri[,4],Mias_halleri[,3]/Mias_halleri[,4],
Piek_halleri[,3]/Piek_halleri[,4],Kowa_halleri[,3]/Kowa_halleri[,4],Zapa_halleri[,3]/Zapa_halleri[,4])

colnames(All_wooutgroup)<-c("Scaffold","Pos","Buko_arenosa_AF","Kato_arenosa_AF","Mias_arenosa_AF",
"Piek_arenosa_AF","Kowa_arenosa_AF","Zapa_arenosa_AF","Buko_halleri_AF","Mias_halleri_AF","Piek_halleri_AF","Kowa_halleri_AF","Zapa_halleri_AF")
All_wooutgroup2<-All_wooutgroup[!(apply(All_wooutgroup[,3:13],1,sum)==0),]
#3770167 of 5180523
All_wooutgroup3<-All_wooutgroup2[complete.cases(All_wooutgroup2),]
#3769786

#HMA2
#scaffold_7	4971705	4976723	7019
HMA2<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=4969705&All3$Pos<=4976723,]
nrow(HMA2)
#269 sites
HMA2_new<-HMA2[apply(HMA2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
nrow(HMA2_new)
#265 SNPs
pdf("Heatmap_SNPs_HMA2.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


#HMA2 plus promoter
#scaffold_7	4969705  4976723	7019
HMA2<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=4969705&All3$Pos<=4976723,]
nrow(HMA2)
#269 sites
HMA2_new<-HMA2[apply(HMA2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
nrow(HMA2_new)
#265 SNPs
pdf("Heatmap_SNPs_HMA2_plus_promoter.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

#HMA3
#scaffold_7	4960535	4963880	3346
HMA3<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=4960535&All3$Pos<=4963880,]
nrow(HMA3)
#289 SNPs
HMA3_new<-HMA3[apply(HMA3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

#HMA4
#scaffold_3	23475941	23485785	9845
HMA4<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=23475941&All3$Pos<=23483785,]
nrow(HMA4)
#63 SNPs
HMA4_new<-HMA4[apply(HMA4[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA4.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA4_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=2.5,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

#MTP1
#scaffold_4	22776960	22778882	1923
MTP1<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=22776960&All3$Pos<=22778882,]
nrow(MTP1)
#95 SNPs
MTP1_new<-MTP1[apply(MTP1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MTP1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MTP1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=2.5,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

#IRT1
#AL7G34350	scaffold_7	10085962	10089516
IRT1<-All[All3$Scaffold=="scaffold_7"&All3$Pos>=10085962&All3$Pos<=10087516,]
nrow(IRT1)
#119 SNPs
IRT1_new<-IRT1[apply(IRT1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_IRT1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(IRT1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

#NRAMP3
#AL4G13010	scaffold_4	1634727	1637073
NRAMP3<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=1634727&All3$Pos<=1637073,]
nrow(NRAMP3)
#127 SNPs
NRAMP3_new<-NRAMP3[apply(NRAMP3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_NRAMP3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(NRAMP3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

require(openxlsx)
gff<-read.delim("D:/Lyrata/Alyrata_384_v2.1.gene.gff3",sep="\t",header=F,comment.char="#")
RBH<-read.xlsx("E:/Lyrata_RBH/Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx",1)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
gff2_Gene <- strsplit(gff$Attributes,";",fixed=TRUE)
gff2_Gene2<-sapply(gff2_Gene,head,1)
gff2_Gene3<-strsplit(gff2_Gene2,"\\.")
gff2_Gene4<-sapply(gff2_Gene3,head,1)
gff2_Gene6<-gsub("ID=","",gff2_Gene4)
gff$Gene<-gff2_Gene6

gff_RBH<-merge(gff[gff$Type=="gene",],RBH,by.x="Gene",by.y="Alyr_ID")
gff_RBH[grep("ZIP6",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]



ZIP6<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=13949254&All3$Pos<=13950804,]
nrow(ZIP6)
#37 SNPs
ZIP6_new<-ZIP6[apply(ZIP6[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIP6.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIP6_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[grep("ZIP11",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
ZIP11<-All3[All3$Scaffold=="scaffold_1"&All3$Pos>=29573151&All3$Pos<=29574444,]
nrow(ZIP11)
#86 SNPs
ZIP11_new<-ZIP11[apply(ZIP11[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIP11.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIP11_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT1G55920",c(1:6,11:13,18)]
SERAT2_1<-All3[All3$Scaffold=="scaffold_1"&All3$Pos>=29579011&All3$Pos<=29580087,]
nrow(SERAT2_1)
#87 SNPs
SERAT2_1_new<-SERAT2_1[apply(SERAT2_1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_SERAT2_1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(SERAT2_1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[gff_RBH$Ath_ID=="AT2G15290",c(1:6,11:13,18)]
PIC1<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=19789770&All3$Pos<=19791689,]
nrow(PIC1)
#76 SNPs
PIC1_new<-PIC1[apply(PIC1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_PIC1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(PIC1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("ZTP29",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
ZTP29<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=9213498&All3$Pos<=9216853,]
nrow(ZTP29)
#256 SNPs
ZTP29_new<-ZTP29[apply(ZTP29[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZTP29.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZTP29_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[gff_RBH$Ath_ID=="AT4G36060",c(1:6,11:13,18)]
bHLH011<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=2048538&All3$Pos<=2050007,]
nrow(bHLH011)
#114 SNPs
bHLH011_new<-bHLH011[apply(bHLH011[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_bHLH011.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(bHLH011_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[grep("SPL7",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
Spl7<-All3[All3$Scaffold=="scaffold_6"&All3$Pos>=7862984&All3$Pos<=7867548,]
nrow(Spl7)
#353 SNPs
Spl7_new<-Spl7[apply(Spl7[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_Spl7.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(Spl7_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[grep("YSL2",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
YSL2<-All3[All3$Scaffold=="scaffold_6"&All3$Pos>=10644250&All3$Pos<=10647305,]
nrow(YSL2)
#229 SNPs
YSL2_new<-YSL2[apply(YSL2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_YSL2.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(YSL2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[gff_RBH$Ath_ID=="AT3G58810",c(1:6,11:13,18)]
MTP3<-All3[All3$Scaffold=="scaffold_5"&All3$Pos>=19047210&All3$Pos<=19049108,]
nrow(MTP3)
#50 SNPs
MTP3_new<-MTP3[apply(MTP3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MTP3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MTP3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[gff_RBH$Gene=="AL7G24770",c(1:6,11:13,18)]
MEB3<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=5991780&All3$Pos<=5996124,]
nrow(MEB3)
#295 SNPs
MEB3_new<-MEB3[apply(MEB3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MEB3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MEB3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("FRO2",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
FRO2<-All3[All3$Scaffold=="scaffold_1"&All3$Pos>=266658&All3$Pos<=270551,]
nrow(FRO2)
#219 SNPs
FRO2_new<-FRO2[apply(FRO2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FRO2.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FRO2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[gff_RBH$Ath_ID=="AT4G40030",c(1:6,11:13,18)]
HTR4<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=876859&All3$Pos<=878213,]
nrow(HTR4)
#60 SNPs
HTR4_new<-HTR4[apply(HTR4[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HTR4.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HTR4_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[gff_RBH$Ath_ID=="AT3G13110",c(1:6,11:13,18)]
SERAT2_2<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=1634727&All3$Pos<=1637073,]
nrow(SERAT2_2)
#127 SNPs
SERAT2_2_new<-SERAT2_2[apply(SERAT2_2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_SERAT2_2.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(SERAT2_2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("HMA1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
#AL7G13570	scaffold_7
HMA1<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=1400419&All3$Pos<=1405046,]
nrow(HMA1)
#365 SNPs
HMA1_new<-HMA1[apply(HMA1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[grep("HMA1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
#AL7G13570	scaffold_7
HMA1<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=1398419&All3$Pos<=1405046,]
nrow(HMA1)
#483 SNPs
HMA1_new<-HMA1[apply(HMA1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA1_withpromoter.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[gff_RBH$Ath_ID=="AT2G23240",c(1:6,11:13,18)]
MT4b<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=1727007&All3$Pos<=1728574,]
nrow(MT4b)
#13 SNPs
MT4b_new<-MT4b[apply(MT4b[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MT4b.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MT4b_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("PAA1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
#AL7G18020	AT4G33520	scaffold_7	3229360	3233354
PAA1<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=3229360&All3$Pos<=3233354,]
nrow(PAA1)
#241 SNPs
PAA1_new<-PAA1[apply(PAA1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_PAA1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(PAA1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("FRO8",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
FRO8<-All3[All3$Scaffold=="scaffold_8"&All3$Pos>=13737963&All3$Pos<=13741359,]
nrow(FRO8)
#360 SNPs
FRO8_new<-FRO8[apply(FRO8[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FRO8.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FRO8_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("ZIP8",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
ZIP8<-All3[All3$Scaffold=="scaffold_2"&All3$Pos>=13440351&All3$Pos<=13441078,]
nrow(ZIP8)
#46 SNPs
ZIP8_new<-ZIP8[apply(ZIP8[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIP8.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIP8_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("CAX2",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
CAX2<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=5374636&All3$Pos<=5378276,]
nrow(CAX2)
#203 SNPs
CAX2_new<-CAX2[apply(CAX2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_CAX2.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(CAX2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("IRT2",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
IRT2<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=10093506&All3$Pos<=10094856,]
nrow(IRT2)
#93 SNPs
IRT2_new<-IRT2[apply(IRT2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_IRT2.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(IRT2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("FER3",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
FER3<-All3[All3$Scaffold=="scaffold_5"&All3$Pos>=17807314&All3$Pos<=17809288,]
nrow(FER3)
#115 SNPs
FER3_new<-FER3[apply(FER3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FER3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FER3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("IRT3",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
IRT3<-All3[All3$Scaffold=="scaffold_2"&All3$Pos>=2822105&All3$Pos<=2824141,]
nrow(IRT3)
#165 SNPs
IRT3_new<-IRT3[apply(IRT3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_IRT3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(IRT3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("FER1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
FER1<-All3[All3$Scaffold=="scaffold_6"&All3$Pos>=187631&All3$Pos<=189383,]
nrow(FER1)
#59 SNPs
FER1_new<-FER1[apply(FER1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FER1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FER1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[gff_RBH$Ath_ID=="AT3G08040",c(1:6,11:13,18)]
FRD3<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=3352004&All3$Pos<=3357699,]
nrow(FRD3)
#179 SNPs
FRD3_new<-FRD3[apply(FRD3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FRD3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FRD3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("ZIF1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
ZIF1<-All3[All3$Scaffold=="scaffold_6"&All3$Pos>=5543272&All3$Pos<=5547484,]
nrow(ZIF1)
#229 SNPs
ZIF1_new<-ZIF1[apply(ZIF1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIF1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIF1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("EIN2",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
EIN2<-All3[All3$Scaffold=="scaffold_6"&All3$Pos>=927173&All3$Pos<=933065,]
nrow(EIN2)
#330 SNPs
EIN2_new<-EIN2[apply(EIN2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_EIN2.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(EIN2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("NRAMP5",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
NRAMP5<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=10603275&All3$Pos<=10605431,]
nrow(NRAMP5)
#194 SNPs
NRAMP5_new<-NRAMP5[apply(NRAMP5[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_NRAMP5.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(NRAMP5_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("ZIFL1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
ZIFL1<-All3[All3$Scaffold=="scaffold_6"&All3$Pos>=5552512&All3$Pos<=5555951,]
nrow(ZIFL1)
#245 SNPs
ZIFL1_new<-ZIFL1[apply(ZIFL1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIFL1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIFL1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("OPT3",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
OPT3<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=12164276&All3$Pos<=12167047,]
nrow(OPT3)
#214 SNPs
OPT3_new<-OPT3[apply(OPT3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_OPT3.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(OPT3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()




gff_RBH[grep("CHR5",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
CHR5<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=18100712&All3$Pos<=18113115,]
nrow(CHR5)
#752 SNPs
CHR5_new<-CHR5[apply(CHR5[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_CHR5.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(CHR5_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

gff_RBH[grep("MOS1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
MOS1<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=7470769&All3$Pos<=7478325,]
nrow(MOS1)
#559 SNPs
MOS1_new<-MOS1[apply(MOS1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MOS1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MOS1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[grep("HMA1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
HMA1_alt<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=1738889&All3$Pos<=1742253,]
nrow(HMA1_alt)
#6 SNPs
HMA1_alt_new<-HMA1_alt[apply(HMA1_alt[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA1_alt.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA1_alt_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[grep("HMA1",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
HMA1_alt<-All_wooutgroup[All_wooutgroup$Scaffold=="scaffold_4"&All_wooutgroup$Pos>=1738889&All_wooutgroup$Pos<=1742253,]
nrow(HMA1_alt)
#8 SNPs
HMA1_alt_new<-HMA1_alt[apply(HMA1_alt[,3:13],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA1_alt_wooutgroup.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA1_alt_new[,rev(c(10,11,9,13,12,5,6,3,4,8,7))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa"))
cols_pops<-c(rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:11)
	{mtext(side=2,at=seq(0,1,length.out=11)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()





#CCA1 and LHY2




gff_RBH[gff_RBH$Ath_ID=="AT2G46830",c(1:6,11:13,18)]
CCA1<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=22785130&All3$Pos<=22788405,]
nrow(CCA1)
#272 SNPs
CCA1_new<-CCA1[apply(CCA1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_CCA1.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(CCA1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT2G46830",c(1:6,11:13,18)]
CCA1<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=22782958&All3$Pos<=22788941,]
nrow(CCA1)
#383 SNPs
CCA1_new<-CCA1[apply(CCA1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_CCA1_2500bpupstream_700bpdownstream_added.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(CCA1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT1G01060",c(1:6,11:13,18)]
LHY<-All3[All3$Scaffold=="scaffold_1"&All3$Pos>=502430&All3$Pos<=504430,]
nrow(LHY)
#90 SNPs
LHY_new<-LHY[apply(LHY[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_LHY_promoter_Twisstwindow.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(LHY_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT2G46790",c(1:6,11:13,18)]
PRR9<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=22767360&All3$Pos<=22769933,]
nrow(PRR9)
#176 SNPs
PRR9_new<-PRR9[apply(PRR9[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_PRR9_Twisstwindow.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(PRR9_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT3G46640",c(1:6,11:13,18)]
PCL1<-All3[All3$Scaffold=="scaffold_5"&All3$Pos>=12733964&All3$Pos<=12736171,]
nrow(PCL1)
#81 SNPs
PCL1_new<-PCL1[apply(PCL1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_PCL1_Twisstwindow.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(PCL1_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[gff_RBH$Ath_ID=="AT2G25930",c(1:6,11:13,18)]
ELF3<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=5420517&All3$Pos<=5424627,]
nrow(ELF3)
#240 SNPs
ELF3_new<-ELF3[apply(ELF3[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ELF3_Twisstwindow.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ELF3_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT4G01120",c(1:6,11:13,18)]
GBF2<-All3[All3$Scaffold=="scaffold_6"&All3$Pos>=24485282&All3$Pos<=24488190,]
nrow(GBF2)
#204 SNPs
GBF2_new<-GBF2[apply(GBF2[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_GBF2_Twisstwindow.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(GBF2_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT2G18790",c(1:6,11:13,18)]
PhyB<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=23698491&All3$Pos<=23703073,]
nrow(PhyB)
#294 SNPs
PhyB_new<-PhyB[apply(PhyB[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_PhyB_Twisstwindow.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(PhyB_new[,rev(c(15,16,14,18,17,10,11,8,9,13,12,7,4,5,6,3))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.09,-0.01,-0.09,0.279, xpd = TRUE,lwd=2)
segments(-0.09,0.3,-0.09,0.69, xpd = TRUE,lwd=2)
segments(-0.09,0.72,-0.09,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

