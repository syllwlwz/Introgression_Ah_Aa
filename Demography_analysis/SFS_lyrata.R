

AFdata_KK<-read.table("HalleriKromKosiGS_new.csv",header=TRUE,sep="\t")
AFdata_KK$AF<-AFdata_KK$AC/AFdata_KK$AN
AFdata_KK$AF.1<-AFdata_KK$AC.1/AFdata_KK$AN.1
names(AFdata_KK)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1","AF","AF.1")
d_KK=AFdata_KK[AFdata_KK$AF<max(AFdata_KK$AF) | AFdata_KK$AF.1<max(AFdata_KK$AF.1),]
d_KK_Krom<-ifelse(d_KK$AF>0.5,1-d_KK$AF,d_KK$AF)
d_KK_Kosi<-ifelse(d_KK$AF.1>0.5,1-d_KK$AF.1,d_KK$AF.1)

AFdata_NP<-read.table("HalleriNossPaisGS_new.csv",header=TRUE,sep="\t")
AFdata_NP$AF<-AFdata_NP$AC/AFdata_NP$AN
AFdata_NP$AF.1<-AFdata_NP$AC.1/AFdata_NP$AN.1
names(AFdata_NP)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1","AF","AF.1")
d_NP=AFdata_NP[AFdata_NP$AF<max(AFdata_NP$AF) | AFdata_NP$AF.1<max(AFdata_NP$AF.1),]
d_NP_Noss<-ifelse(d_NP$AF>0.5,1-d_NP$AF,d_NP$AF)
d_NP_Pais<-ifelse(d_NP$AF.1>0.5,1-d_NP$AF.1,d_NP$AF.1)

AFdata_LB<-read.table("HalleriLangBestGS_new.csv",header=TRUE,sep="\t")
AFdata_LB$AF<-AFdata_LB$AC/AFdata_LB$AN
AFdata_LB$AF.1<-AFdata_LB$AC.1/AFdata_LB$AN.1
names(AFdata_LB)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1","AF","AF.1")
d_LB=AFdata_LB[AFdata_LB$AF<max(AFdata_LB$AF) | AFdata_LB$AF.1<max(AFdata_LB$AF.1),]
d_LB_Lang<-ifelse(d_LB$AF>0.5,1-d_LB$AF,d_LB$AF)
d_LB_Best<-ifelse(d_LB$AF.1>0.5,1-d_LB$AF.1,d_LB$AF.1)

pdf("SFS_KKNPLB.pdf", width=18, height=6, paper="special",pointsize=15)
par(mfrow=c(1,3))
par(mar=c(1,1,1,1)+0.1)
par(oma=c(4,5,0,0))
par(mgp=c(5,0.2,0))
options(scipen=10)

maxKr<-hist(d_KK_Krom,breaks=10,plot=F)
maxKo<-hist(d_KK_Kosi,breaks=10,plot=F)
hist(d_KK_Krom,col=rgb(t(col2rgb("red")),alpha=150,maxColorValue=255),xlab=expression(bold(Allele~frequency~Krom-Kosi)),ylab=expression(bold(Frequency)),main="",breaks=10,cex=1,cex.lab=1,cex.axis=1,mgp=c(2,1,0),ylim=c(0,1500000),axes=F)
hist(d_KK_Kosi,col=rgb(t(col2rgb("grey20")),alpha=150,maxColorValue=255),add=TRUE,breaks=10)
legend("topright",legend=c(expression(Krom),expression(Kosi)),cex=2.5,pch=15,col=rgb(t(col2rgb(c("red","grey20"))),alpha=150,maxColorValue=255),pt.cex=3,bty="n",inset=0.05)
axis(1,line=0,cex.axis=1.5,cex.lab=2,padj=0.75)
axis(2,at=c(0,400000,800000,1200000),line=0,cex.axis=1.5,cex.lab=2,mgp=c(5,0.5,0),labels=formatC(c(0,400000,800000,1200000)/1000,format="d",big.mark=","))
box(lwd=1.2)
maxNo<-hist(d_NP_Noss,breaks=10,plot=F)
maxPa<-hist(d_NP_Pais,breaks=10,plot=F)
hist(d_NP_Noss,col=rgb(t(col2rgb("red")),alpha=150,maxColorValue=255),xlab=expression(bold(Allele~frequency~Noss-Pais)),ylab=expression(bold(Frequency)),main="",breaks=10,cex=1,cex.lab=1,cex.axis=1,mgp=c(2,1,0),ylim=c(0,1500000),axes=F)
hist(d_NP_Pais,col=rgb(t(col2rgb("grey20")),alpha=150,maxColorValue=255),add=TRUE,breaks=10)
legend("topright",legend=c(expression(Noss),expression(Pais)),cex=2.5,pch=15,col=rgb(t(col2rgb(c("red","grey20"))),alpha=150,maxColorValue=255),pt.cex=3,bty="n",inset=0.05)
axis(1,line=0,cex.axis=1.5,cex.lab=2,padj=0.75)
axis(2,at=c(0,400000,800000,1200000),line=0,cex.axis=1.5,cex.lab=2,mgp=c(5,0.5,0),labels=formatC(c(0,400000,800000,1200000)/1000,format="d",big.mark=","))
box(lwd=1.2)

maxLa<-hist(d_LB_Lang,breaks=10,plot=F)
maxBe<-hist(d_LB_Best,breaks=10,plot=F)
hist(d_LB_Lang,col=rgb(t(col2rgb("red")),alpha=150,maxColorValue=255),xlab=expression(bold(Allele~frequency~Lang-Best)),ylab=expression(bold(Frequency)),main="",breaks=10,cex=1,cex.lab=1,cex.axis=1,mgp=c(2,1,0),ylim=c(0,1500000),axes=F)
hist(d_LB_Best,col=rgb(t(col2rgb("grey20")),alpha=150,maxColorValue=255),add=TRUE,breaks=10)
legend("topright",legend=c(expression(Lang),expression(Best)),cex=2.5,pch=15,col=rgb(t(col2rgb(c("red","grey20"))),alpha=150,maxColorValue=255),pt.cex=3,bty="n",inset=0.05)
axis(1,line=0,cex.axis=1.5,cex.lab=2,padj=0.75)
axis(2,at=c(0,400000,800000,1200000),line=0,cex.axis=1.5,cex.lab=2,mgp=c(5,0.5,0),labels=formatC(c(0,400000,800000,1200000)/1000,format="d",big.mark=","))

mtext(side=2,line=2,"Density [mille]",outer=T,cex=1.8)
mtext(side=1,line=2,"Frequency",outer=T,cex=1.8)
box(lwd=1.2)
dev.off()

require(DescTools)
pdf("SFS_KKNPLB_orderered.pdf", width=12, height=4, paper="special",pointsize=12)
par(mfrow=c(1,3))
par(mar=c(1,1,1,1)+0.1)
par(oma=c(4,5,0,0))
par(mgp=c(5,0.2,0))
options(scipen=10)
par(lwd=2)

barplot(t(cbind(Freq(d_KK_Kosi,breaks=10)$perc,Freq(d_KK_Krom,breaks=10)$perc)),main="",cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("black","red"),col="white",space=c(0.1,1),ylim=c(0,0.5),las=1)
legend("topright",legend=c("Krom","Kosi"),pch=15,fill="white",cex=1.5,bty="n",border=c("red","black"),inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)

barplot(t(cbind(Freq(d_NP_Pais,breaks=10)$perc,Freq(d_NP_Noss,breaks=10)$perc)),main="",cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("black","red"),col="white",space=c(0.1,1),ylim=c(0,0.5),las=1)
legend("topright",legend=c("Noss","Pais"),pch=15,fill="white",cex=1.5,bty="n",border=c("red","black"),inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)

barplot(t(cbind(Freq(d_LB_Best,breaks=10)$perc,Freq(d_LB_Lang,breaks=10)$perc)),main="",cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("black","red"),col="white",space=c(0.1,1),ylim=c(0,0.5),las=1)
legend("topright",legend=c("Lang","Best"),pch=15,fill="white",cex=1.5,bty="n",border=c("red","black"),inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)

mtext(side=1,line=2,"Allele frequency",outer=T,cex=1.2)
mtext(side=2,line=2,"Proportion of SNPs",outer=T,cex=1.2)

dev.off()


