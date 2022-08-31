library("SNPRelate")
vcf.fn<-"All5_introgression_4dgsites.LD_Pruned.vcf"
snpgdsVCF2GDS(vcf.fn,"All5_introgression_4dgsites.gds",method="biallelic.only",ignore.chr.prefix=c("scaffold_"))
genofile<-snpgdsOpen("All5_introgression_4dgsites.gds")
#showfile.gds(closeall=TRUE) if left open
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
#Number of samples: 186
#import 23301 variants


Introgression_pca<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples)
#identifies monomorphic SNPs probably because no MAF and AN very high, so beyond decimal places of R
#Principal Component Analysis (PCA) on genotypes:
    # of samples: 186
    # of SNPs: 23301
#    using 8 threads
    # of principal components: 32
#PCA:    the sum of all selected genotypes (0,1,2) = 6691681


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
,rep("green",26)#Buko
,rep("blue",19)#Kato
,rep("orange",19)#Kowa
,rep("purple",34)#Mias
,rep("lightblue",34)#Piek
,rep("yellow",26)#Zapa
)


pdf("PCA_Introgression_woutgroup.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
dev.off()


Introgression_pca_sub<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(1:28)])
#    # of samples: 158
    # of SNPs: 23,301

require(openxlsx)
write.xlsx(cbind(Introgression_pca_sub$sample.id,Introgression_pca_sub$ eigenvect[,1]),"PCA_All5_PC1.xlsx")

ev <- Introgression_pca_sub$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#3 PCs

colour_pops_sub<-c(rep("orange",26)#Buko
,rep("magenta",19)#Kato
,rep("black",19)#Kowa
,rep("red",34)#Mias
,rep("purple",34)#Piek
,rep("grey",26)#Zapa
)




pdf("PCA_Introgression2.pdf",width=10,height=10,paper="special",pointsize=16)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"))
legend("topright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"))
dev.off()

pdf("PCA_Introgression.pdf",width=8,height=6,paper="special",pointsize=12)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], pch=17,col=colour_pops_sub,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_sub)[c(3,6,2,1,5,4)],legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa")[c(3,6,2,1,5,4)],cex=1.3,border=c("black"),inset=c(-0.16,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),cex=1.5)
dev.off()

pdf("PCA_Introgression_wooutgroup_LDpruned_MAF_0_05_PC1_3.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,3], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC3 (",round(Introgression_pca_sub$varprop[3]*100,0),"%)"),cex=1.5)
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,3], labels=Introgression_pca_sub$sample, cex= 0.4)
dev.off()

Genotyping_ploidy4<-read.table("../Total_samples_introgression_4.list",sep=",",header=F)
Genotyping_ploidy2<-read.table("../Total_samples_introgression_2.list",sep=",",header=F)
GT_ploidy4<-cbind(Genotyping_ploidy4,rep(4,length(Genotyping_ploidy4)))
colnames(GT_ploidy4)<-c("Sample","Ploidy")
GT_ploidy2<-cbind(Genotyping_ploidy2,rep(2,length(Genotyping_ploidy2)))
colnames(GT_ploidy2)<-c("Sample","Ploidy")
Missing<-data.frame(c("Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43","Piek_a20_1","Piek_a20_2","Piek_h14","Piek_h6"),c(2,2,2,2,4,4,4,4))
colnames(Missing)<-c("Sample","Ploidy")


Genotyping_ploidy<-rbind(GT_ploidy4,GT_ploidy2,Missing)
Ploidy<-merge(as.data.frame(Introgression_pca_sub$sample),Genotyping_ploidy,by.x="Introgression_pca_sub$sample",by.y="Sample",all.x=T)
colnames(Ploidy)<-c("Sample","Ploidy")


pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2.pdf",width=15,height=13,paper="special",pointsize=30)
par(mar=c(5,5,1,8))
par(mgp=c(2,0.65,0))
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_sub$sample,Ploidy$Sample)])],col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),las=1,cex=1.5)
legend("bottomright",fill=unique(colour_pops_sub)[c(3,6,2,1,5,4)],legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa")[c(3,6,2,1,5,4)],cex=1.3,border=c("black"),inset=c(-0.3,0),xpd=T,bty="n")
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n",inset=c(-0.45,0),xpd=T)
dev.off()

pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2_PC1_3.pdf",width=15,height=15,paper="special",pointsize=25)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,3],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_sub$sample,Ploidy$Sample)])],col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC3 (",round(Introgression_pca_sub$varprop[3]*100,0),"%)"))
legend("topright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),bty="n",cex=1.3)
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,3], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,inset=c(-0.18,0),border=c("black"),bty="n")
dev.off()

pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2_PC2_3_names.pdf",width=15,height=13,paper="special",pointsize=30)
par(mar=c(5,5,1,8))
par(mgp=c(2,0.65,0))
plot(Introgression_pca_sub$eigenvect[,2],Introgression_pca_sub$eigenvect[,3],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_sub$sample,Ploidy$Sample)])],col=colour_pops_sub,xlab=paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),ylab=paste0("PC3 (",round(Introgression_pca_sub$varprop[3]*100,0),"%)"),las=1,cex=1.5)
legend("bottomright",fill=unique(colour_pops_sub)[c(3,6,2,1,5,4)],legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa")[c(3,6,2,1,5,4)],cex=1.3,border=c("black"),inset=c(-0.3,0),xpd=T,bty="n")
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n",inset=c(-0.45,0),xpd=T)
dev.off()

pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2_PC2_4_ns.pdf",width=15,height=15,paper="special",pointsize=25)
plot(Introgression_pca_sub$eigenvect[,2],Introgression_pca_sub$eigenvect[,4],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_sub$sample,Ploidy$Sample)])],col=colour_pops_sub,xlab=paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),ylab=paste0("PC4 (",round(Introgression_pca_sub$varprop[4]*100,0),"%)"))
legend("bottomright",fill=unique(colour_pops_sub),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),bty="n",cex=1.3)
#text(Introgression_pca_sub$eigenvect[,2],Introgression_pca_sub$eigenvect[,4], labels=Introgression_pca_sub$sample, cex= 0.7)
#text(Introgression_pca_sub$eigenvect[,2][c(22,24,38,68,70,87,129,157,167)],Introgression_pca_sub$eigenvect[,3][c(22,24,38,68,70,87,129,157,167)], labels=Introgression_pca_sub$sample[c(22,24,38,68,70,87,129,157,167)], cex= 0.7)
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n")
dev.off()




pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2_horizleg.pdf",width=16,height=15,paper="special",pointsize=30)
par(mar=c(5,5,4,1))
par(mgp=c(4,0.65,0))
par(lwd=3)
plot(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_sub$sample,Ploidy$Sample)])],col=colour_pops_sub,xlab=paste("PC1 (",round(Introgression_pca_sub$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),las=1,cex=1.3,cex.lab=1.8,cex.axis=1.5)
legend("topleft",fill=unique(colour_pops_sub)[c(3,6,2,1,5,4)],legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa")[c(3,6,2,1,5,4)],cex=1,border=c("black"),inset=c(-0.1,-0.12),xpd=T,bty="n",horiz=T)
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("topleft",pch=c(17,19),legend=c("Diploid","Tetraploid"),cex=1,border=c("black"),bty="n",inset=c(-0.09,-0.2),xpd=T,horiz=T)
dev.off()

pdf("PCA_Introgression2_LDpruned_MAF_0_05_version2_PC2_3_names_horizleg.pdf",width=16,height=15,paper="special",pointsize=30)
par(mar=c(5,5,4,1))
par(mgp=c(4,0.65,0))
par(lwd=3)
plot(Introgression_pca_sub$eigenvect[,2],Introgression_pca_sub$eigenvect[,3],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_sub$sample,Ploidy$Sample)])],col=colour_pops_sub,xlab=paste0("PC2 (",round(Introgression_pca_sub$varprop[2]*100,0),"%)"),ylab=paste0("PC3 (",round(Introgression_pca_sub$varprop[3]*100,0),"%)"),las=1,cex=1.3,cex.lab=1.8,cex.axis=1.5)
legend("topleft",fill=unique(colour_pops_sub)[c(3,6,2,1,5,4)],legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa")[c(3,6,2,1,5,4)],cex=1,border=c("black"),inset=c(-0.1,-0.12),xpd=T,bty="n",horiz=T)
#text(Introgression_pca_sub$eigenvect[,1],Introgression_pca_sub$eigenvect[,2], labels=Introgression_pca_sub$sample, cex= 0.7)
legend("topleft",pch=c(17,19),legend=c("Diploid","Tetraploid"),cex=1,border=c("black"),bty="n",inset=c(-0.09,-0.2),xpd=T,horiz=T)
dev.off()




Introgression_pca_noAth<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[-c(21:28)])
# of SNPs: 22,495

ev <- Introgression_pca_noAth$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#4 PCs

colour_pops_noAth<-c(rep("red",20)#Outgroup
,rep("green",26)#Buko
,rep("blue",19)#Kato
,rep("orange",19)#Kowa
,rep("purple",34)#Mias
,rep("lightblue",34)#Piek
,rep("yellow",26)#Zapa
)

Genotyping_ploidy<-rbind(GT_ploidy4,GT_ploidy2,Missing)
Ploidy<-merge(as.data.frame(Introgression_pca_noAth$sample),Genotyping_ploidy,by.x="Introgression_pca_noAth$sample",by.y="Sample",all.x=T)
colnames(Ploidy)<-c("Sample","Ploidy")
Ploidy$Ploidy[1:20]<-2



pdf("PCA_Introgression_noAth.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_noAth$sample,Ploidy$Sample)])],col=colour_pops_noAth,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"),cex=1.5)
text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.4)
dev.off()

pdf("PCA_Introgression_noAth_PC3_4.pdf",width=14,height=12,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_noAth$eigenvect[,3],Introgression_pca_noAth$eigenvect[,4],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_noAth$sample,Ploidy$Sample)])],col=colour_pops_noAth,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC3 (",round(Introgression_pca_noAth$varprop[3]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC4 (",round(Introgression_pca_noAth$varprop[4]*100,0),"%)"),cex=1.5)
#text(Introgression_pca_noAth$eigenvect[,3],Introgression_pca_noAth$eigenvect[,4],labels=Introgression_pca_noAth$sample, cex= 0.4)
text(-0.05,0.145,labels="A. croatica",font=3)
text(0.2,0.13,labels="A. lyrata",font=3)
text(0.29,-0.2,labels="A. pedemontana",font=3)
text(0.28,-0.27,labels="A. cebennensis",font=3)
legend("bottomleft",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.01,0),xpd=T,bty="n")
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n")
dev.off()


pdf("PCA_Introgression_noAth_PC1_2_labelled.pdf",width=14,height=12,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_noAth$sample,Ploidy$Sample)])],col=colour_pops_noAth,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomleft",fill=unique(colour_pops_noAth),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.01,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"),cex=1.5)
#text(Introgression_pca_noAth$eigenvect[,1],Introgression_pca_noAth$eigenvect[,2], labels=Introgression_pca_noAth$sample, cex= 0.4)
text(0.05,-0.3,labels="A. croatica",font=3)
text(0.05,0.09,labels="A. lyrata",font=3)
text(0.017,-0.2,labels="A. pedemontana",font=3)
text(0.02,-0.22,labels="A. cebennensis",font=3)
legend("topleft",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n")
dev.off()


pdf("PCA_Introgression_woAth_LDpruned_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
#text(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], labels=Introgression_pca$sample, cex= 0.4)
dev.off()









#############
#only lyrata#
#############

Introgression_pca_lyr<-snpgdsPCA(genofile,num.thread=8,algorithm="exact",autosome.only=FALSE,sample.id=samples[-c(19:28,1:9)])
## of SNPs: 21,931

ev <- Introgression_pca_lyr$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
#3 PCs

colour_pops_lyr<-c(rep("red",9)#Outgroup
,rep("green",26)#Buko
,rep("blue",19)#Kato
,rep("orange",19)#Kowa
,rep("purple",34)#Mias
,rep("lightblue",34)#Piek
,rep("yellow",26)#Zapa
)

Genotyping_ploidy<-rbind(GT_ploidy4,GT_ploidy2,Missing)
Ploidy<-merge(as.data.frame(Introgression_pca_lyr$sample),Genotyping_ploidy,by.x="Introgression_pca_lyr$sample",by.y="Sample",all.x=T)
colnames(Ploidy)<-c("Sample","Ploidy")
Ploidy$Ploidy[1:9]<-2

pdf("PCA_Introgression_lyr_PC3_4.pdf",width=14,height=12,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_lyr$eigenvect[,3],Introgression_pca_lyr$eigenvect[,4],pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_lyr$sample,Ploidy$Sample)])],col=colour_pops_lyr,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_lyr),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.4,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC3 (",round(Introgression_pca_lyr$varprop[3]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC4 (",round(Introgression_pca_lyr$varprop[4]*100,0),"%)"),cex=1.5)
#text(Introgression_pca_noAth$eigenvect[,3],Introgression_pca_lyr$eigenvect[,4],labels=Introgression_pca_lyr$sample, cex= 0.4)
legend("bottomleft",fill=unique(colour_pops_noAth),legend=c("A. lyrata","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.01,0),xpd=T,bty="n")
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n")
dev.off()


pdf("PCA_Introgression_lyr_PC1_2_labelled.pdf",width=14,height=12,paper="special",pointsize=16)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca_lyr$eigenvect[,1],Introgression_pca_lyr$eigenvect[,2], pch=c(17,19)[as.factor(Ploidy$Ploidy[match(Introgression_pca_lyr$sample,Ploidy$Sample)])],col=colour_pops_lyr,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops_lyr),legend=c("A. lyrata","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca_noAth$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca_noAth$varprop[2]*100,0),"%)"),cex=1.5)
#text(Introgression_pca_lyr$eigenvect[,1],Introgression_pca_lyr$eigenvect[,2], labels=Introgression_pca_lyr$sample, cex= 0.4)
legend("topright",pch=c(17,19),legend=c("Diploid","Tetraploid"),title="Genotyped as",cex=1.3,border=c("black"),bty="n")
dev.off()











pdf("PCA_Introgression_woAth_LDpruned_MAF_0_05.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(Introgression_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(Introgression_pca$varprop[2]*100,0),"%)"),cex=1.5)
#text(Introgression_pca$eigenvect[,1],Introgression_pca$eigenvect[,2], labels=Introgression_pca$sample, cex= 0.4)
dev.off()








#IBD

genofile<-snpgdsOpen("All5_introgression_4dgsites.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

#Estimating IBD Using PLINK method of moments (MoM)

# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile,maf=0,missing.rate=1,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples)
# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

#Identity-By-State Analysis
#For n study individuals, snpgdsIBS() can be used to create a n×n matrix of genome-wide average IBS pairwise identities:
#without outgroup
#Identity-By-State Analysis
#For n study individuals, snpgdsIBS() can be used to create a n×n matrix of genome-wide average IBS pairwise identities:

ibs2 <- snpgdsIBS(genofile,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples[-c(1:28)],maf=0,missing.rate=1)
ibs_df2<-data.frame(ibs2$ibs)
names(ibs_df2)<-ibs2$sample.id
colnames(ibs_df2)<-ibs2$sample.id
rownames(ibs_df2)<-ibs2$sample.id

Genotyping_ploidy4<-read.table("../Total_samples_introgression_4.list",sep=",",header=F)
Genotyping_ploidy2<-read.table("../Total_samples_introgression_2.list",sep=",",header=F)
GT_ploidy4<-cbind(Genotyping_ploidy4,rep(4,length(Genotyping_ploidy4)))
colnames(GT_ploidy4)<-c("Sample","Ploidy")
GT_ploidy2<-cbind(Genotyping_ploidy2,rep(2,length(Genotyping_ploidy2)))
colnames(GT_ploidy2)<-c("Sample","Ploidy")
Missing<-data.frame(c("Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43","Piek_a20_1","Piek_a20_2","Piek_h14","Piek_h6"),c(2,2,2,2,4,4,4,4))
colnames(Missing)<-c("Sample","Ploidy")

Genotyping_ploidy<-rbind(GT_ploidy4,GT_ploidy2,Missing)
Ploidy<-merge(as.data.frame(rownames(ibs_df2)),Genotyping_ploidy,by.x="rownames(ibs_df2)",by.y="Sample",all.x=T)
colnames(Ploidy)<-c("Sample","Ploidy")
Ploidy$Pop<-substr(Ploidy$Sample,1,4)
Ploidy$Pop[Ploidy$Pop=="Zako"]<-"Zapa"
summary1<-aggregate(ibs_df2,list(Ploidy$Pop,Ploidy$Ploidy),mean)

summary2<-aggregate(t(summary1[,-c(1:2)]),list(Ploidy$Pop,Ploidy$Ploidy),mean)
colnames(summary2)[3:14]<-paste(summary1$Group.1,summary1$Group.2,sep="_")
rownames(summary2)<-c(paste(summary1$Group.1,summary1$Group.2,sep="_"))
write.xlsx(summary2[,-c(1:2)],"IBS_populations_species_introgression.xlsx",row.names=T,overwrite=T)










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
require(openxlsx)
write.xlsx(sub_ibs3,"IBS_above_0_9_Introgression_wooutgroup_LD_pruned_MAF_0_05.xlsx",overwrite=T)

write.xlsx(sub_ibs2,"IBS_above_0_95_Introgression_wooutgroup_LD_pruned_MAF_0_05.xlsx",overwrite=T)
write.xlsx(ibs_df2,"IBS_Introgression_wooutgroup_LD_pruned_MAF_0_05.xlsx",row.names=T,overwrite=T)

#with outgroup, MAF 0.05
ibs4 <- snpgdsIBS(genofile,num.thread=2,autosome.only=FALSE,remove.monosnp=F,sample.id=samples,maf=0,missing.rate=1)
ibs_df4<-data.frame(ibs4$ibs)
names(ibs_df4)<-ibs4$sample.id
colnames(ibs_df4)<-ibs4$sample.id
rownames(ibs_df4)<-ibs4$sample.id

pdf("Hist_IBS_Introgression_woutgroup_LD_pruned_MAF_0_05.pdf",width=10,height=8,paper="special",pointsize=16)
hist(ibs4$ibs)
dev.off()

write.xlsx(ibs_df4,"IBS_Introgression_woutgroup_LD_pruned_MAF_0_05.xlsx",row.names=T,overwrite=T)



Genotyping_ploidy4<-read.table("../Total_samples_introgression_4.list",sep=",",header=F)
Genotyping_ploidy2<-read.table("../Total_samples_introgression_2.list",sep=",",header=F)
GT_ploidy4<-cbind(Genotyping_ploidy4,rep(4,length(Genotyping_ploidy4)))
colnames(GT_ploidy4)<-c("Sample","Ploidy")
GT_ploidy2<-cbind(Genotyping_ploidy2,rep(2,length(Genotyping_ploidy2)))
colnames(GT_ploidy2)<-c("Sample","Ploidy")
Missing<-data.frame(c("Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43","Piek_a20_1","Piek_a20_2","Piek_h14","Piek_h6"),c(2,2,2,2,4,4,4,4))
colnames(Missing)<-c("Sample","Ploidy")

Genotyping_ploidy<-rbind(GT_ploidy4,GT_ploidy2,Missing)
Ploidy<-merge(as.data.frame(rownames(ibs_df4)),Genotyping_ploidy,by.x="rownames(ibs_df4)",by.y="Sample",all.x=T)
colnames(Ploidy)<-c("Sample","Ploidy")
Ploidy$Pop<-substr(Ploidy$Sample,1,4)
Ploidy$Pop[Ploidy$Pop=="Zako"]<-"Zapa"
Ploidy$Ploidy[1:28]<-2

summary1<-aggregate(ibs_df4,list(Ploidy$Pop,Ploidy$Ploidy),mean)

summary2<-aggregate(t(summary1[,-c(1:2)]),list(Ploidy$Pop,Ploidy$Ploidy),mean)
colnames(summary2)[3:19]<-paste(summary1$Group.1,summary1$Group.2,sep="_")
rownames(summary2)<-c(paste(summary1$Group.1,summary1$Group.2,sep="_"))
require(openxlsx)
write.xlsx(summary2[,-c(1:2)],"IBS_populations_species_introgression_woutgroup.xlsx",row.names=T,overwrite=T)






