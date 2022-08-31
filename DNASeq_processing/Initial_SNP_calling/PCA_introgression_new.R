library("SNPRelate")
vcf.fn<-"All_introgression_4dgsites.vcf"
snpgdsVCF2GDS(vcf.fn,"All_introgression_4dgsites.gds",method="biallelic.only",ignore.chr.prefix=c("scaffold_"))
genofile<-snpgdsOpen("All_introgression_4dgsites.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

Introgression_pca<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-60])
#identifies monomorphic SNPs probably because no MAF and AN very high, so beyond decimal places of R


jpeg("Eigenvalues_Kaiser-Guttman-Test.jpeg", width=26, height=18, units="cm", res=1000)
ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
barplot (ev, main ="Eigenvalues", col="black", las=2)
abline (h=mean(ev,na.rm=T), col="red")
dev.off()

pdf("Eigenvalues_Kaiser-Guttman-Test.pdf",width=8,height=6,paper="special",pointsize=16)
ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
barplot (ev, main ="Eigenvalues", col="black", las=2)
abline (h=mean(ev,na.rm=T), col="red")
dev.off()



ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#5 PCs

colour_pops<-c(rep("red",28)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)


pdf("PCA_Introgression_woutgroup.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
dev.off()


jpeg("PCA_Introgression.jpeg", width=26, height=18, units="cm", res=1000)

plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=17,col=colour_pops,xlab=paste("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"))
legend("bottomright",fill=unique(colour_pops),legend=c("Buko","Kowa","Outgroup","Kato","Mias","Piek","Zapa"))
text(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2],labels=Introgression_pca$sample, cex= 0.7)


dev.off()


plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,3], pch=17,col=colour_pops,xlab=paste("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),ylab=paste("PC3 (",round(Introgression_pca$varprop[3]*100,0),"%)"))
legend("bottomright",fill=unique(colour_pops),legend=c("Buko","Kowa","Outgroup","Kato","Mias","Piek","Zapa"))
text(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,3], labels=Introgression_pca$sample, cex= 0.7)


plot(Introgression_pca$eigenvect[,2],Introgression_pca$eigenvect[,3], pch=17,col=colour_pops,xlab=paste("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),ylab=paste("PC3 (",round(Introgression_pca$varprop[3]*100,0),"%)"))
legend("bottomright",fill=unique(colour_pops),legend=c("Buko","Kowa","Outgroup","Kato","Mias","Piek","Zapa"))
text(Introgression_pca$eigenvect[,2],Introgression_pca$eigenvect[,3], labels=Introgression_pca$sample, cex= 0.7)



Introgression_pca_sub<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(1:28,60)])

ev <- Introgression_pca_sub$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#4 PCs

colour_pops_sub<-c(rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)


pdf("PCA_Introgression2.pdf",width=20,height=20,paper="special",pointsize=16)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"))
legend("topleft",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"))
dev.off()

pdf("PCA_Introgression.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
dev.off()









Introgression_pca_noAth<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[-c(21:28,60)])

ev <- Introgression_pca_noAth$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#6 PCs

colour_pops_noAth<-c(rep("red",20)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)



jpeg("PCA_Introgression_noAth.jpeg", width=26, height=18, units="cm", res=1000)
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], pch=17,col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"))
legend("bottomleft",fill=unique(colour_pops_noAth),legend=c("Buko","Kowa","Kato","Mias","Piek","Zapa"))
text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.7)
dev.off()

pdf("PCA_Introgression_noAth.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], pch=17,col=colour_pops_noAth,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.4)
dev.off()



#for Vero: Mias Zapa
Introgression_pca_MZ<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[c(97:132,170:205)])

ev <- Introgression_pca_MZ$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#4 PCs

colour_pops_MZ<-c(rep("red",36)#Mias
,rep("black",36)#Zapa
)

pdf("PCA_Introgression_MZ.pdf",width=10,height=8,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_MZ$eigenvect[,1],Introgression_pca_MZ$eigenvect[,2], pch=17,col=colour_pops_MZ,xlab="",ylab="",cex=1.25,cex.lab=1.5,las=1,cex.axis=1.25)
legend("bottomright",fill=unique(colour_pops_MZ),legend=c("Mias","Zapa"),cex=1.3,border=c("black"),inset=c(-0.2,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_MZ$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_MZ$varprop[2]*100,0),"%)"),cex=1.5)
#text(Introgression_pca_MZ$eigenvect[,1],Introgression_pca_MZ$eigenvect[,2], labels=Introgression_pca_MZ$sample, cex=0.4)
dev.off()

#############################
#with MAF 0.05#
#############################
library("SNPRelate")
vcf.fn<-"All_introgression_4dgsites.vcf"
snpgdsVCF2GDS(vcf.fn,"All_introgression_4dgsites.gds",method="biallelic.only",ignore.chr.prefix=c("scaffold_"))
genofile<-snpgdsOpen("All_introgression_4dgsites.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

Introgression_pca<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-60],maf=0.05)
#Excluding 923,586 SNPs (monomorphic: FALSE, MAF: 0.05, missing rate: NaN)
#Working space: 204 samples, 269,622 SNPs

ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#5 PCs

colour_pops<-c(rep("red",28)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)

pdf("PCA_Introgression_woutgroup_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
dev.off()

Introgression_pca_sub<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(1:28,60)],maf=0.05)

ev <- Introgression_pca_sub$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#4 PCs

colour_pops_sub<-c(rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)

pdf("PCA_Introgression2_MAF_0_05.pdf",width=20,height=20,paper="special",pointsize=16)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"))
legend("topleft",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"))
dev.off()

pdf("PCA_Introgression_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
dev.off()

pdf("PCA_Introgression_samplenames_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
dev.off()

pdf("PCA_Introgression_samplenames_subset_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1][c(22,24,38,129,68,70,157,167)],Introgression_pca_sub$eigenvect[,2][c(22,24,38,129,68,70,157,167)], labels=Introgression_pca_sub$sample[c(22,24,38,129,68,70,157,167)], cex= 0.7)
dev.off()

Introgression_pca_noAth<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[-c(21:28,60)],maf=0.05)
#Excluding 934,455 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: NaN)
#Working space: 196 samples, 258,753 SNPs

ev <- Introgression_pca_noAth$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#3 PCs

colour_pops_noAth<-c(rep("red",20)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)

pdf("PCA_Introgression_noAth_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], pch=17,col=colour_pops_noAth,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.4)
dev.off()



#############################
#LD pruned#
#############################
library("SNPRelate")
vcf.fn<-"All_introgression_4dgsites.LD_Pruned.vcf"
snpgdsVCF2GDS(vcf.fn,"All_introgression_4dgsites.LD_Pruned.gds",method="biallelic.only",ignore.chr.prefix=c("scaffold_"))
genofile<-snpgdsOpen("All_introgression_4dgsites.LD_Pruned.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

Introgression_pca<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-60])
#identifies monomorphic SNPs probably because no MAF and AN very high, so beyond decimal places of R


jpeg("Eigenvalues_Kaiser-Guttman-Test.jpeg", width=26, height=18, units="cm", res=1000)
ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
barplot (ev, main ="Eigenvalues", col="black", las=2)
abline (h=mean(ev,na.rm=T), col="red")
dev.off()

pdf("Eigenvalues_Kaiser-Guttman-Test.pdf",width=8,height=6,paper="special",pointsize=16)
ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
barplot (ev, main ="Eigenvalues", col="black", las=2)
abline (h=mean(ev,na.rm=T), col="red")
dev.off()



ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#5 PCs

colour_pops<-c(rep("red",28)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)


pdf("PCA_Introgression_woutgroup_LDpruned.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
dev.off()

Introgression_pca_sub<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(1:28,60,97)])

ev <- Introgression_pca_sub$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#4 PCs

colour_pops_sub<-c(rep("orange",27)#Buko
,rep("magenta",19)#Kato
,rep("black",21)#Kowa
,rep("red",36)#Mias
,rep("purple",37)#Piek
,rep("grey",36)#Zapa
)


pdf("PCA_Introgression2_LDpruned.pdf",width=20,height=20,paper="special",pointsize=16)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"))
legend("topright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"))
dev.off()

pdf("PCA_Introgression_LDpruned.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
dev.off()


pdf("PCA_Introgression_LDpruned_samplenames_duplicates_version2.pdf",width=15,height=15,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1][c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175)],Introgression_pca_sub$eigenvect[,2][c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175)], labels=Introgression_pca_sub$sample[c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175)], cex= 0.7)
points(Introgression_pca_sub$eigenvect[c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175),1],Introgression_pca_sub$eigenvect[c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175),2], pch=17,col="red",xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
dev.off()

pdf("PCA_Introgression_LDpruned_samplenames.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
dev.off()


pdf("PCA_Introgression_LDpruned_samplenames_subset.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1][c(22,24,38,129,68,70,157,167)],Introgression_pca_sub$eigenvect[,2][c(22,24,38,129,68,70,157,167)], labels=Introgression_pca_sub$sample[c(22,24,38,129,68,70,157,167)], cex= 0.7)
dev.off()


Introgression_pca_noAth<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[-c(21:28,60)])

ev <- Introgression_pca_noAth$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#7 PCs

colour_pops_noAth<-c(rep("red",20)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)



jpeg("PCA_Introgression_noAth_LDpruned.jpeg", width=26, height=18, units="cm", res=1000)
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], pch=17,col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"))
legend("bottomleft",fill=unique(colour_pops_noAth),legend=c("Buko","Kowa","Kato","Mias","Piek","Zapa"))
text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.7)
dev.off()

pdf("PCA_Introgression_noAth_LDpruned.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], pch=17,col=colour_pops_noAth,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.4)
dev.off()

#for Vero: Mias Zapa
Introgression_pca_MZ<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[c(97:132,170:205)])

ev <- Introgression_pca_MZ$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#4 PCs

colour_pops_MZ<-c(rep("red",36)#Mias
,rep("black",36)#Zapa
)

pdf("PCA_Introgression_MZ_LDpruned.pdf",width=10,height=8,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_MZ$eigenvect[,1],Introgression_pca_MZ$eigenvect[,2], pch=17,col=colour_pops_MZ,xlab="",ylab="",cex=1.25,cex.lab=1.5,las=1,cex.axis=1.25)
legend("bottomright",fill=unique(colour_pops_MZ),legend=c("Mias","Zapa"),cex=1.3,border=c("black"),inset=c(-0.2,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_MZ$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_MZ$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_MZ$eigenvect[,1],Introgression_pca_MZ$eigenvect[,2], labels=Introgression_pca_MZ$sample, cex=0.4)
dev.off()


#IBD

genofile<-snpgdsOpen("All_introgression_4dgsites.LD_Pruned.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

#Estimating IBD Using PLINK method of moments (MoM)

# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile,maf=0,missing.rate=1,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-60])
# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

#Estimating IBD Using Maximum Likelihood Estimation (MLE)

# Estimate IBD coefficients
ibd <- snpgdsIBDMLE(genofile,maf=0, missing.rate=1,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-60])

#Identity-By-State Analysis
#For n study individuals, snpgdsIBS() can be used to create a n×n matrix of genome-wide average IBS pairwise identities:

ibs <- snpgdsIBS(genofile,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-60])
ibs_df<-data.frame(ibs$ibs)
names(ibs_df)<-ibs$sample.id
colnames(ibs_df)<-ibs$sample.id
rownames(ibs_df)<-ibs$sample.id

indx <- ibs_df>0.95
rn <- rownames(ibs_df)[row(ibs_df)*indx]
cn <-  colnames(ibs_df)[col(ibs_df)*indx]
val <- ibs_df[indx]
sub_ibs<-data.frame(rn, cn, val)

#without outgroup
#Identity-By-State Analysis
#For n study individuals, snpgdsIBS() can be used to create a n×n matrix of genome-wide average IBS pairwise identities:

ibs2 <- snpgdsIBS(genofile,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(1:28,60)],maf=0.05,missing.rate=1)
ibs_df2<-data.frame(ibs2$ibs)
names(ibs_df2)<-ibs2$sample.id
colnames(ibs_df2)<-ibs2$sample.id
rownames(ibs_df2)<-ibs2$sample.id

indx2 <- ibs_df2>0.95
rn2 <- rownames(ibs_df2)[row(ibs_df2)*indx2]
cn2 <-  colnames(ibs_df2)[col(ibs_df2)*indx2]
val2 <- ibs_df2[indx2]
sub_ibs2<-data.frame(rn2, cn2, val2)

pdf("Hist_IBS_Introgression_wooutgroup_LD_pruned_MAF_0_05.pdf",width=10,height=8,paper="special",pointsize=16)
hist(ibs2$ibs)
dev.off()

indx3 <- ibs_df2>0.9
rn3 <- rownames(ibs_df2)[row(ibs_df2)*indx3]
cn3 <-  colnames(ibs_df2)[col(ibs_df2)*indx3]
val3 <- ibs_df2[indx3]
sub_ibs3<-data.frame(rn3, cn3, val3)
write.xlsx(sub_ibs3,"IBS_above_0_9_Introgression_wooutgroup_LD_pruned_MAF_0_05.xlsx")

write.xlsx(sub_ibs2,"IBS_above_0_95_Introgression_wooutgroup_LD_pruned_MAF_0_05.xlsx")
write.xlsx(ibs_df2,"IBS_Introgression_wooutgroup_LD_pruned_MAF_0_05.xlsx",row.names=T)

write.xlsx(ibs_df,"IBS_Introgression_LD_pruned.xlsx",row.names=T)

pdf("Hist_IBS_Introgression_wooutgroup_LD_pruned_woutgroup.pdf",width=10,height=8,paper="special",pointsize=16)
hist(ibs$ibs)
dev.off()

#with outgroup, MAF 0.05
ibs4 <- snpgdsIBS(genofile,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(60)],maf=0.05,missing.rate=1)
#Excluding 25,919 SNPs (monomorphic: FALSE, MAF: 0.05, missing rate: 1)
#Working space: 204 samples, 7,215 SNPs
ibs_df4<-data.frame(ibs4$ibs)
names(ibs_df4)<-ibs4$sample.id
colnames(ibs_df4)<-ibs4$sample.id
rownames(ibs_df4)<-ibs4$sample.id

pdf("Hist_IBS_Introgression_woutgroup_LD_pruned_MAF_0_05.pdf",width=10,height=8,paper="special",pointsize=16)
hist(ibs4$ibs)
dev.off()

write.xlsx(ibs_df4,"IBS_Introgression_woutgroup_LD_pruned_MAF_0_05.xlsx",row.names=T)



#############################
#LD pruned with MAF 0.05#
#############################
library("SNPRelate")
vcf.fn<-"All_introgression_4dgsites.LD_Pruned.vcf"
snpgdsVCF2GDS(vcf.fn,"All_introgression_4dgsites.LD_Pruned2.gds",method="biallelic.only",ignore.chr.prefix=c("scaffold_"))
genofile<-snpgdsOpen("All_introgression_4dgsites.LD_Pruned2.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

Introgression_pca<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-60],maf=0.05)
#identifies monomorphic SNPs probably because no MAF and AN very high, so beyond decimal places of R
#Excluding 25,919 SNPs (monomorphic: FALSE, MAF: 0.05, missing rate: NaN)
#Working space: 204 samples, 7,215 SNPs

ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#5 PCs

colour_pops<-c(rep("red",28)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)

Genotyping_ploidy<-read.table("nQuire_HC_final_old.list",sep=",",header=F)

pdf("PCA_Introgression_woutgroup_LDpruned_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=c(17,19)[as.factor(Genotyping_ploidy$V2)],col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
#text(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], labels=Introgression_pca$sample, cex= 0.7)
legend("bottomright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),inset=c(-0.45,0.65),xpd=T,bty="n")
dev.off()

Introgression_pca_sub<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(1:28,60,97)],maf=0.05)

ev <- Introgression_pca_sub$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#3 PCs

colour_pops_sub<-c(rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)

pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2.pdf",width=20,height=20,paper="special",pointsize=25)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"))
legend("topright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"))
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("bottomright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n")
dev.off()

pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2_PC1_3.pdf",width=20,height=20,paper="special",pointsize=25)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,3],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC3 (",round(Introgression_pca_sub$varprop[3]*100,0),"%)"))
legend("topleft",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"))
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,3], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,inset=c(0.45,0),border=c("black"),bty="n")
dev.off()

pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2_PC2_3_names.pdf",width=20,height=20,paper="special",pointsize=25)
plot(Introgression_pca_sub$eigenvect[,2],Introgression_pca_sub$eigenvect[,3],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),ylab=paste("PC3 (",round(Introgression_pca_sub$varprop[3]*100,0),"%)"))
legend("topleft",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"))
#text(Introgression_pca_sub$eigenvect[,2],Introgression_pca_sub$eigenvect[,3], labels=Introgression_pca_sub$sample, cex= 0.7)
#text(Introgression_pca_sub$eigenvect[,2][c(22,24,38,68,70,87,129,157,167)],Introgression_pca_sub$eigenvect[,3][c(22,24,38,68,70,87,129,157,167)], labels=Introgression_pca_sub$sample[c(22,24,38,68,70,87,129,157,167)], cex= 0.7)
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n")
dev.off()







pdf("PCA_Introgression_LDpruned_MAF_0_05_version2.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("bottomright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),inset=c(-0.45,0.5),xpd=T,bty="n")
dev.off()

pdf("PCA_Introgression_LDpruned_samplenames_MAF_0_05_version2.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("bottomright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),inset=c(-0.45,0.5),xpd=T,bty="n")
dev.off()

pdf("PCA_Introgression_LDpruned_samplenames_subset_MAF_0_05_version2.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1][c(22,24,38,129,68,70,87,157,167)],Introgression_pca_sub$eigenvect[,2][c(22,24,38,129,68,70,87,157,167)], labels=Introgression_pca_sub$sample[c(22,24,38,129,68,70,157,167)], cex= 0.7)
text(Introgression_pca_sub$eigenvect[,1][c(69)],Introgression_pca_sub$eigenvect[,2][c(69)], labels=Introgression_pca_sub$sample[c(69)], cex= 0.7)
points(Introgression_pca_sub$eigenvect[164,1],Introgression_pca_sub$eigenvect[164,2], pch=17,col="red",xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),inset=c(-0.45,0.5),xpd=T,bty="n")
dev.off()

pdf("PCA_Introgression_LDpruned_samplenames_duplicates_MAF_0_05_version2.pdf",width=15,height=15,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1][c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175)],Introgression_pca_sub$eigenvect[,2][c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175)], labels=Introgression_pca_sub$sample[c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175)], cex= 0.7)
points(Introgression_pca_sub$eigenvect[c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175),1],Introgression_pca_sub$eigenvect[c(1,6,60:63,126,127,137,138,141,142,150,152,155:159,160:161,164:169,173,175),2], pch=17,col="red",xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),inset=c(-0.45,0.5),xpd=T,bty="n")
dev.off()

pdf("PCA_Introgression_LDpruned_samplenames_subset_MAF_0_05_Piek_h6clone_version2.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(1:28)])],col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_sub$eigenvect[,1][c(140)],Introgression_pca_sub$eigenvect[,2][c(140)], labels=Introgression_pca_sub$sample[c(140)], cex= 0.7)
legend("bottomright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),inset=c(-0.45,0.5),xpd=T,bty="n")
dev.off()


Introgression_pca_noAth<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[-c(21:28,60)],maf=0.05)

ev <- Introgression_pca_noAth$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#7 PCs

colour_pops_noAth<-c(rep("red",20)#Outgroup
,rep("green",27)#Buko
,rep("blue",19)#Kato
,rep("orange",21)#Kowa
,rep("purple",36)#Mias
,rep("lightblue",37)#Piek
,rep("yellow",36)#Zapa
)

pdf("PCA_Introgression_noAth_LDpruned_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2],pch=c(17,19)[as.factor(Genotyping_ploidy$V2[-c(21:28)])],col=colour_pops_noAth,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.4)
dev.off()


#############################
#LD pruned with MAF 0.05 2#
#############################
library("SNPRelate")
vcf.fn<-"All_introgression_4dgsites.LD_Pruned.MAF.vcf"
snpgdsVCF2GDS(vcf.fn,"All_introgression_4dgsites.LD_Pruned.MAF.gds",method="biallelic.only",ignore.chr.prefix=c("scaffold_"))
genofile<-snpgdsOpen("All_introgression_4dgsites.LD_Pruned.MAF.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
#Number of samples: 176
#Parsing "All_introgression_4dgsites.LD_Pruned.MAF.vcf" ...
#        import 10525 variants.

Introgression_pca<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-68],maf=0)

ev <- Introgression_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#3 PCs

require(StAMPP)
require(vcfR)
require(adegenet)

vcf<-read.vcfR("All_introgression_4dgsites.LD_Pruned.MAF.vcf", verbose = FALSE)
##2.2 Convert vcf into genlight object
#leads to NA for tetraploids
Intro_GL <- vcfR2genlight.tetra(vcf)
Ploidy<-ploidy(Intro_GL)[-68]

colour_pops<-c(rep("orange",27)#Buko
,rep("magenta",19)#Kato
,rep("black",21)#Kowa
,rep("red",35)#Mias
,rep("purple",37)#Piek
,rep("grey",36)#Zapa
)


pdf("PCA_Introgression_wooutgroup_LDpruned_MAF_0_05_paper2.pdf",width=12,height=8,paper="special",pointsize=16)
par(mar=c(5,5,5,1))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=c(17,19)[as.factor(Ploidy)],col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("topleft",fill=unique(colour_pops),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.1,-0.2),xpd=T,bty="n",horiz=T)
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
#text(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], labels=Introgression_pca$sample, cex= 0.4)
dev.off()

pdf("PCA_Introgression_wooutgroup_LDpruned_MAF_0_05_prefil_PC1_3.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,3], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC3 (",round(Introgression_pca$varprop[3]*100,0),"%)"),cex=1.5)
text(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,3], labels=Introgression_pca$sample, cex= 0.4)
dev.off()


