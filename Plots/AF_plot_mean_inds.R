options(java.parameters = "-Xmx12g")
require(openxlsx)

All<-read.table("All5_repolarized_AFpersample.vcf",sep="\t",header=F)
colnames(All)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","A_ce_1","A_ce_2","A_ce_3","A_cr_1","A_cr_2","A_cr_3","A_cr_4","A_cr_5","A_cr_6","A_ly_1","A_ly_2","A_ly_3","A_ly_4","A_ly_5","A_ly_6","A_ly_7",
"A_ly_ICE1","A_ly_Plech","A_pe_1","A_pe_2","A_th_1","A_th_2","A_th_3","A_th_4","A_th_5","A_th_6","A_th_7","A_th_8","Buko_a21","Buko_a22","Buko_a28","Buko_a29","Buko_a30","Buko_a31","Buko_a34","Buko_a34b","Buko_h04","Buko_h12",
"Buko_h1_7","Buko_h21","Buko_h22","Buko_h23","Buko_h24","Buko_h25","Buko_h26","Buko_h27","Buko_h28","Buko_h29","Buko_h30","Buko_h31","Buko_h32","Buko_h33","Buko_h35","Buko_h36","Kato_a21","Kato_a22","Kato_a23","Kato_a24",
"Kato_a27","Kato_a28","Kato_a29","Kato_a30","Kato_a33","Kato_a35","Kato_h09","Kato_h21","Kato_h22","Kato_h23","Kato_h24","Kato_h26","Kato_h27","Kato_h33","Kato_h35","Kowa001a04","Kowa001a05","Kowa001a06","Kowa001a07",
"Kowa001a08","Kowa001a09","Kowa001a11","Kowa001a12","Kowa002h13","Kowa002h14","Kowa002h17","Kowa002h18","Kowa002h21","Kowa_h02comb","Kowa_h04comb","Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43","Mias001a18","Mias002h19","Mias002h20",
"Mias003a09","Mias003a10","Mias003a11","Mias003a13","Mias003a15","Mias003a16","Mias004h02","Mias004h06","Mias004h07","Mias004h08","Mias009h03","Mias009h05","Mias_19a03x21a01_d","Mias_23a03x22a03_k","Mias_a12xa05_l","Mias_a42",
"Mias_a43","Mias_a44","Mias_a45","Mias_a46","Mias_a47","Mias_a48","Mias_a50","Mias_a54","Mias_a55","Mias_h09","Mias_h41","Mias_h43","Mias_h44","Mias_h45","Mias_h46","Piek_a1","Piek_a10","Piek_a12","Piek_a13","Piek_a14",
"Piek_a15","Piek_a16","Piek_a17","Piek_a18","Piek_a19","Piek_a2","Piek_a20_1","Piek_a20_2","Piek_a22","Piek_a4","Piek_a5","Piek_a7","Piek_a8","Piek_a9","Piek_h01","Piek_h03","Piek_h04","Piek_h05comb","Piek_h06comb","Piek_h07",
"Piek_h10","Piek_h11","Piek_h12","Piek_h13","Piek_h14","Piek_h14_2","Piek_h15","Piek_h2","Piek_h6","Zako002h02","Zako002h03","Zako002h04","Zako002h05","Zako002h06","Zako002h07","Zako002h08","Zako_h01comb","Zako_h03",
"Zako_h09comb","Zako_h10","Zako_h12","Zapa002a01","Zapa002a02","Zapa002a03","Zapa002a04","Zapa002a05","Zapa002a07","Zapa004a11","Zapa008a19","Zapa_11a03x12a01_l","Zapa_a09xa03_b","Zapa_h01comb","Zapa_h07","Zapa_h11comb",
"Zapa_h12")

for (i in 10:195)
	{All[,i]<-as.numeric(All[,i])
	}


Pop_species<-c(rep("Format",9),paste(substr(colnames(All)[10:195],1,4),c(rep(2,28),rep(4,8),rep(2,18)
,rep(4,10),2,rep(4,16),rep(2,11),4,2,2,rep(4,6),rep(2,6),rep(4,11),2,4,rep(2,6),rep(4,19),rep(2,10),4,rep(2,16),rep(4,10),rep(2,4)),sep="_"))
Pop_species[Pop_species=="Zako_2"]<-"Zapa_2"

Pop_species_levels<-levels(as.factor(Pop_species))[!(levels(as.factor(Pop_species))=="Kato_2"|levels(as.factor(Pop_species))=="Format")]

All2<-data.frame(All[,1:9])
count=10
for (j in Pop_species_levels)
	{All2[,count]<-apply(All[,Pop_species==j],1,function(x) {median(x[which(!is.na(x))])})
	count=count+1
	}

All2[,24]<-apply(All[,Pop_species=="Zapa_2"],1,function(x) {median(x[which(!is.na(x))])})
colnames(All2)[24]<-"Zapa_2"

colnames(All2)[10:26]<-Pop_species_levels


All3<-All2[!(apply(All2[,15:25],1,sum)==0),]
All3<-All3[complete.cases(All3[,15:25]),]
#1435940
require(RColorBrewer)

colnames(All3)[1:2]<-c("Scaffold","Pos")
#HMA2
#scaffold_7	4969705	4976723	7019
HMA2<-All2[All3$Scaffold=="scaffold_7"&All3$Pos>=4971705&All3$Pos<=4976723,]
nrow(HMA2)
#417 sites
HMA2_new<-HMA2[apply(HMA2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
nrow(HMA2_new)
#118 SNPs
pdf("Heatmap_SNPs_HMA2_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#scaffold_7	4971705  4976723	7019
HMA2<-All3[All3$Scaffold=="scaffold_7"&All3$Pos>=4969705&All3$Pos<=4976723,]
nrow(HMA2)
#139 sites
HMA2_new<-HMA2[apply(HMA2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
nrow(HMA2_new)
#126 SNPs
pdf("Heatmap_SNPs_HMA2_plus_promoter_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#157 SNPs
HMA3_new<-HMA3[apply(HMA3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.1,-0.01,-0.1,0.279, xpd = TRUE,lwd=2)
segments(-0.1,0.3,-0.1,0.69, xpd = TRUE,lwd=2)
segments(-0.1,0.72,-0.1,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()

#HMA4
#scaffold_3	23475941	23485785	9845
HMA4<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=23475941&All3$Pos<=23483785,]
nrow(HMA4)
#28 SNPs
HMA4_new<-HMA4[apply(HMA4[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA4_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA4_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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



HMA4<-All3[All3$Scaffold=="scaffold_3"&All3$Pos>=23475941&All3$Pos<=23485785,]
nrow(HMA4)
#31 SNPs
HMA4_new<-HMA4[apply(HMA4[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA4_medianind_withpromoter.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA4_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#54 SNPs
MTP1_new<-MTP1[apply(MTP1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MTP1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MTP1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#120 SNPs
IRT1_new<-IRT1[apply(IRT1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_IRT1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(IRT1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#63 SNPs
NRAMP3_new<-NRAMP3[apply(NRAMP3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_NRAMP3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(NRAMP3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.11,-0.01,-0.11,0.279, xpd = TRUE,lwd=2)
segments(-0.11,0.3,-0.11,0.69, xpd = TRUE,lwd=2)
segments(-0.11,0.72,-0.11,1.01, xpd = TRUE,lwd=2)

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
#15 SNPs
ZIP6_new<-ZIP6[apply(ZIP6[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIP6_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIP6_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.15,-0.01,-0.15,0.279, xpd = TRUE,lwd=2)
segments(-0.15,0.3,-0.15,0.69, xpd = TRUE,lwd=2)
segments(-0.15,0.72,-0.15,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()


gff_RBH[grep("ZIP11",gff_RBH$Araport11_short_name),c(1:6,11:13,18)]
ZIP11<-All3[All3$Scaffold=="scaffold_1"&All3$Pos>=29573151&All3$Pos<=29574444,]
nrow(ZIP11)
#40 SNPs
ZIP11_new<-ZIP11[apply(ZIP11[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIP11_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIP11_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

pops<-rev(c("Mias","Piek","Buko","Zapa","Kowa","Mias","Piek","Buko","Kato","Zapa","Kowa","A_cr","A_ly","A_ce","A_pe","A_th"))
cols_pops<-c(rep("black",5),rep("black",2),rep("red",4),rep("black",2),rep("red",4))

for (i in 1:16)
	{mtext(side=2,at=seq(0,1,length.out=16)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.13,-0.01,-0.13,0.279, xpd = TRUE,lwd=2)
segments(-0.13,0.3,-0.13,0.69, xpd = TRUE,lwd=2)
segments(-0.13,0.72,-0.13,1.01, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Arabidopsis\noutgroup",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

box()
dev.off()



gff_RBH[gff_RBH$Ath_ID=="AT1G55920",c(1:6,11:13,18)]
SERAT2_1<-All3[All3$Scaffold=="scaffold_1"&All3$Pos>=29579011&All3$Pos<=29580087,]
nrow(SERAT2_1)
#37 SNPs
SERAT2_1_new<-SERAT2_1[apply(SERAT2_1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_SERAT2_1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(SERAT2_1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#41 SNPs
PIC1_new<-PIC1[apply(PIC1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_PIC1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(PIC1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#107 SNPs
ZTP29_new<-ZTP29[apply(ZTP29[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZTP29_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZTP29_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#50 SNPs
bHLH011_new<-bHLH011[apply(bHLH011[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_bHLH011_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(bHLH011_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#164 SNPs
Spl7_new<-Spl7[apply(Spl7[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_Spl7_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(Spl7_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#93 SNPs
YSL2_new<-YSL2[apply(YSL2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_YSL2_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(YSL2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#23 SNPs
MTP3_new<-MTP3[apply(MTP3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MTP3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MTP3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#144 SNPs
MEB3_new<-MEB3[apply(MEB3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MEB3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MEB3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#101 SNPs
FRO2_new<-FRO2[apply(FRO2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FRO2_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FRO2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#27 SNPs
HTR4_new<-HTR4[apply(HTR4[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HTR4_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HTR4_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#63 SNPs
SERAT2_2_new<-SERAT2_2[apply(SERAT2_2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_SERAT2_2_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(SERAT2_2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#126 SNPs
HMA1_new<-HMA1[apply(HMA1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#175 SNPs
HMA1_new<-HMA1[apply(HMA1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_HMA1_withpromoter_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(HMA1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#4 SNPs
MT4b_new<-MT4b[apply(MT4b[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_MT4b_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(MT4b_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#124 SNPs
PAA1_new<-PAA1[apply(PAA1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_PAA1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(PAA1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#116 SNPs
FRO8_new<-FRO8[apply(FRO8[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FRO8_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FRO8_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#16 SNPs
ZIP8_new<-ZIP8[apply(ZIP8[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIP8_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIP8_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#75 SNPs
CAX2_new<-CAX2[apply(CAX2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_CAX2_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(CAX2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#54 SNPs
IRT2_new<-IRT2[apply(IRT2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_IRT2_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(IRT2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#30 SNPs
FER3_new<-FER3[apply(FER3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FER3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FER3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#62 SNPs
IRT3_new<-IRT3[apply(IRT3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_IRT3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(IRT3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#27 SNPs
FER1_new<-FER1[apply(FER1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FER1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FER1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#70 SNPs
FRD3_new<-FRD3[apply(FRD3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_FRD3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(FRD3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#104 SNPs
ZIF1_new<-ZIF1[apply(ZIF1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIF1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIF1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#127 SNPs
EIN2_new<-EIN2[apply(EIN2[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_EIN2_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(EIN2_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#74 SNPs
NRAMP5_new<-NRAMP5[apply(NRAMP5[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_NRAMP5_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(NRAMP5_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#120 SNPs
ZIFL1_new<-ZIFL1[apply(ZIFL1[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_ZIFL1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(ZIFL1_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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
#99 SNPs
OPT3_new<-OPT3[apply(OPT3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_OPT3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(OPT3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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


	

NRAMP3<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=1634727&All3$Pos<=1637073,]
nrow(NRAMP3)
#63 SNPs
NRAMP3_new<-NRAMP3[apply(NRAMP3[,15:25],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_NRAMP3_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
par(mar=c(2,10,1,1))
image(as.matrix(NRAMP3_new[,rev(c(20,22,15,24,18,21,23,16,17,25,19,11,12,10,13,14))]),col=colorRampPalette(brewer.pal(9,"Greys"))(11),ylab="",xlab="",xaxt="none",yaxt="none")

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




#CCA1 and LHY2




gff_RBH[gff_RBH$Ath_ID=="AT2G46830",c(1:6,11:13,18)]
CCA1<-All3[All3$Scaffold=="scaffold_4"&All3$Pos>=22785130&All3$Pos<=22788405,]
nrow(CCA1)
#272 SNPs
CCA1_new<-CCA1[apply(CCA1[,8:18],1,function(x) (sum(x)<11)&(sum(x)>0)),]
pdf("Heatmap_SNPs_CCA1_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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
pdf("Heatmap_SNPs_CCA1_2500bpupstream_700bpdownstream_added_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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
pdf("Heatmap_SNPs_LHY_promoter_Twisstwindow_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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
pdf("Heatmap_SNPs_PRR9_Twisstwindow_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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
pdf("Heatmap_SNPs_PCL1_Twisstwindow_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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
pdf("Heatmap_SNPs_ELF3_Twisstwindow_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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
pdf("Heatmap_SNPs_GBF2_Twisstwindow_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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
pdf("Heatmap_SNPs_PhyB_Twisstwindow_medianind.pdf",width=15,height=20,paper="special",pointsize=20)
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

