#Stampp
require(StAMPP)
require(vcfR)
require(adegenet)
# e.g. fst <- stamppFst(geno)

#from Patrick Monnahan: 1.6.21
# convert to genlight 
vcf<-read.vcfR("All5_introgression_4dgsites.LD_Pruned.vcf", verbose = FALSE)
##2.2 Convert vcf into genlight object
#leads to NA for tetraploids
Intro_GL <- vcfR2genlight.tetra(vcf)
locNames(Intro_GL) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(Intro_GL)<-gsub("Zako","Zapa",substr(indNames(Intro_GL),1,4))           # add pop names: here pop names are first 3 chars of ind name

#check
Intro_GL
#186 genotypes,  23,301 binary SNPs, size: 4.5 Mb
# 140445 (3.24 %) missing data

indNames(Intro_GL)
ploidy(Intro_GL)
pop(Intro_GL)

##2.3 Convert genlight to allele frequencies for StAMPP
Intro_StAMPP<-stamppConvert(Intro_GL, type = "genlight")
Pop_species<-paste0(Intro_StAMPP$pop.names,Intro_StAMPP$ploidy)
for (i in levels(as.factor(Pop_species)))
	{assign(paste0("Mean_frequs_",i),apply(Intro_StAMPP[Pop_species==i,6:23306],2,function(x) sum(x,na.rm=T)/length(x[!is.na(x)])))}

All<-data.frame(as.data.frame(vcf@fix[,1:2]),Mean_frequs_A_ce2,Mean_frequs_A_cr2,Mean_frequs_A_ly2,Mean_frequs_A_pe2,Mean_frequs_A_th2,
Mean_frequs_Buko2,Mean_frequs_Buko4,Mean_frequs_Kato2,Mean_frequs_Kato4,Mean_frequs_Kowa2,Mean_frequs_Kowa4,
Mean_frequs_Mias2,Mean_frequs_Mias4,Mean_frequs_Piek2,Mean_frequs_Piek4,Mean_frequs_Zapa2,Mean_frequs_Zapa4)
colnames(All)<-c("Scaffold","Pos",levels(as.factor(Pop_species)))
All$Pos<-as.numeric(All$Pos)
for (j in 3:19)
	{All[,j]<-as.numeric(All[,j])}

#Repolarize
repol<-read.table("Introgression_07repolarized.lookupKey.minInd_4.txt",sep="\t",header=F)
remove<-read.table("Repol_exclude.txt",sep="\t",header=F)

All_repol1<-All[!(paste0(All$Scaffold,All$Pos)%in%paste0(remove$V1,remove$V2)),]
#21413

All_repol<-All_repol1
All_repol[which(paste0(All_repol$Scaffold,All_repol$Pos)%in%paste0(repol$V1,repol$V2)),3:19]<-1-All_repol[which(paste0(All_repol$Scaffold,All_repol$Pos)%in%paste0(repol$V1,repol$V2)),3:19]

All_repol_arenosa1<-All_repol[,c(1:2,9,11,13,15,17,19)]
All_repol_halleri1<-All_repol[,c(1:2,8,12,14,16,18)]

Sum_arenosa<-apply(All_repol_arenosa1[,3:8],1,sum)
Sum_halleri<-apply(All_repol_halleri1[,3:7],1,sum)
All_repol_arenosa<-All_repol_arenosa1[!(Sum_arenosa==6|Sum_arenosa==0),]
#20280
All_repol_halleri<-All_repol_halleri1[!(Sum_halleri==5|Sum_halleri==0),]
#7230

Sum_all_repol<-apply(All_repol[,c(8:9,11:19)],1,sum)
All_repol_sub<-All_repol[!(Sum_all_repol==11|Sum_all_repol==0),c(8:9,11:19)]
All_frequ<-t(cbind(Freq(All_repol_sub[,10],breaks=10)$perc,Freq(All_repol_sub[,4],breaks=10)$perc,Freq(All_repol_sub[,1],breaks=10)$perc,
Freq(All_repol_sub[,8],breaks=10)$perc,Freq(All_repol_sub[,6],breaks=10)$perc,Freq(All_repol_sub[,11],breaks=10)$perc,
Freq(All_repol_sub[,5],breaks=10)$perc,Freq(All_repol_sub[,3],breaks=10)$perc,Freq(All_repol_sub[,2],breaks=10)$perc,
Freq(All_repol_sub[,9],breaks=10)$perc,Freq(All_repol_sub[,7],breaks=10)$perc))


colour_pops_arenosa<-c("grey",#Zapa,
"black",#Kowa
"magenta",#Kato
"orange",#Buko
"purple",#Piek
"red"#Mias
)

colour_pops_halleri<-c("grey",#Zapa
"black",#Kowa
"orange",#Buko
"purple",#Piek
"red"#Mias
)

arenosa<-t(cbind(Freq(All_repol_arenosa[,8],breaks=10)$perc,Freq(All_repol_arenosa[,5],breaks=10)$perc,Freq(All_repol_arenosa[,4],breaks=10)$perc,
Freq(All_repol_arenosa[,3],breaks=10)$perc,Freq(All_repol_arenosa[,7],breaks=10)$perc,Freq(All_repol_arenosa[,6],breaks=10)$perc))

halleri<-t(cbind(Freq(All_repol_halleri[,7],breaks=10)$perc,Freq(All_repol_halleri[,4],breaks=10)$perc,Freq(All_repol_halleri[,3],breaks=10)$perc,
Freq(All_repol_halleri[,6],breaks=10)$perc,Freq(All_repol_halleri[,5],breaks=10)$perc))




require(DescTools)
pdf("SFS_introgression_orderered_subset_within_species_repolarized_paper.pdf", width=15, height=12, paper="special",pointsize=25)
par(mfrow=c(1,2))
par(mar=c(2.5,1,1,1)+0.1)
par(oma=c(3,3,1,1))
par(mgp=c(5,0.2,0))
options(scipen=10)
par(lwd=2)

barplot(arenosa,main=expression(italic(A.~arenosa)),cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=colour_pops_arenosa,col=colour_pops_arenosa,space=c(0.1,1),ylim=c(0,1),las=1)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=2,at=seq(4.25,71.75,7.5),labels=F)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)
text(x=seq(4.25,71.75,7.5),y=-0.03,labels=round(seq(0,1,(1/9)),2),xpd=T,srt=45,cex=1.2,adj=1)

barplot(halleri,main=expression(italic(A.~halleri)),cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=colour_pops_halleri,col=colour_pops_halleri,space=c(0.1,1),ylim=c(0,1),las=1)
legend("topright",legend=c("Zapa","Kowa","Kato","Buko","Piek","Mias"),pch=15,fill=colour_pops_arenosa,cex=1.2,bty="n",border=colour_pops_arenosa,inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=2,at=seq(3.7,61.3,6.4),labels=F)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)
text(x=seq(3.7,61.3,6.4),y=-0.03,labels=round(seq(0,1,(1/9)),2),xpd=T,srt=45,cex=1.2,adj=1)

mtext(side=1,line=0,"Allele frequency",outer=T,cex=1.2)
mtext(side=2,line=1,"Proportion of SNPs",outer=T,cex=1.2)

dev.off()


pdf("SFS_introgression_orderered_plot_within_species_repolarized_paper.pdf", width=15, height=12, paper="special",pointsize=25)
par(mfrow=c(1,2))
par(mar=c(2.5,1,1,1)+0.1)
par(oma=c(3,3,1,1))
par(mgp=c(5,0.2,0))
options(scipen=10)
par(lwd=2)

barplot(All_frequ[6:11,],main=expression(italic(A.~arenosa)),cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=colour_pops_arenosa,col=colour_pops_arenosa,space=c(0.1,1),ylim=c(0,1),las=1)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=2,at=seq(4.25,71.75,7.5),labels=F)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)
text(x=seq(4.25,71.75,7.5),y=-0.03,labels=round(seq(0,1,(1/9)),2),xpd=T,srt=45,cex=1.2,adj=1)

barplot(All_frequ[1:5,],main=expression(italic(A.~halleri)),cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=colour_pops_halleri,col=colour_pops_halleri,space=c(0.1,1),ylim=c(0,1),las=1)
legend("topright",legend=c("Zapa","Kowa","Kato","Buko","Piek","Mias"),pch=15,fill=colour_pops_arenosa,cex=1.2,bty="n",border=colour_pops_arenosa,inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=2,at=seq(3.7,61.3,6.4),labels=F)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)
text(x=seq(3.7,61.3,6.4),y=-0.03,labels=round(seq(0,1,(1/9)),2),xpd=T,srt=45,cex=1.2,adj=1)

mtext(side=1,line=0,"Allele frequency",outer=T,cex=1.2)
mtext(side=2,line=1,"Proportion of SNPs",outer=T,cex=1.2)

dev.off()

pdf("SFS_introgression_orderered_subset_all_repolarized_paper.pdf", width=15, height=12, paper="special",pointsize=25)
par(mar=c(2.5,1,1,1)+0.1)
par(oma=c(3,3,1,1))
par(mgp=c(5,0.2,0))
options(scipen=10)
par(lwd=3)

barplot(All_frequ,cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("darkgreen","darkblue"),col=c(colour_pops_halleri,colour_pops_arenosa),space=c(0.1,1),ylim=c(0,1),las=1)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=2,at=seq(7,124,13),labels=F)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)
text(x=seq(7,124,13),y=-0.03,labels=round(seq(0,1,(1/9)),2),xpd=T,srt=45,cex=1.2,adj=1)

legend("topright",legend=c("Zapa","Kowa","Kato","Buko","Piek","Mias"),pch=15,fill=colour_pops_arenosa,cex=1.2,bty="n",border=colour_pops_arenosa,inset=0.05,col=NA)

mtext(side=1,line=0,"Allele frequency",outer=T,cex=1.2)
mtext(side=2,line=1,"Proportion of SNPs",outer=T,cex=1.2)

dev.off()

All_repol_unfil<-All_repol[,c(8:9,11:19)]
All_frequ_unfil<-t(cbind(Freq(All_repol_unfil[,10],breaks=10)$perc,Freq(All_repol_unfil[,4],breaks=10)$perc,Freq(All_repol_unfil[,1],breaks=10)$perc,
Freq(All_repol_unfil[,8],breaks=10)$perc,Freq(All_repol_unfil[,6],breaks=10)$perc,Freq(All_repol_unfil[,11],breaks=10)$perc,
Freq(All_repol_unfil[,5],breaks=10)$perc,Freq(All_repol_unfil[,3],breaks=10)$perc,Freq(All_repol_unfil[,2],breaks=10)$perc,
Freq(All_repol_unfil[,9],breaks=10)$perc,Freq(All_repol_unfil[,7],breaks=10)$perc))


pdf("SFS_introgression_orderered_subset_all_repolarized_paper_unfiltered.pdf", width=15, height=12, paper="special",pointsize=25)
par(mar=c(2.5,1,1,1)+0.1)
par(oma=c(3,3,1,1))
par(mgp=c(5,0.2,0))
options(scipen=10)
par(lwd=3)

barplot(All_frequ_unfil,cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("darkgreen","darkblue"),col=c(colour_pops_halleri,colour_pops_arenosa),space=c(0.1,1),ylim=c(0,1),las=1)
axis(1,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=2,at=seq(7,124,13),labels=F)
axis(2,line=0,cex.axis=1.2,cex.lab=1.3,mgp=c(5,0.5,0),las=1)
box(lwd=1)
text(x=seq(7,124,13),y=-0.03,labels=round(seq(0,1,(1/9)),2),xpd=T,srt=45,cex=1.2,adj=1)

legend("topright",legend=c("Zapa","Kowa","Kato","Buko","Piek","Mias"),pch=15,fill=colour_pops_arenosa,cex=1.2,bty="n",border=colour_pops_arenosa,inset=0.05,col=NA)

mtext(side=1,line=0,"Allele frequency",outer=T,cex=1.2)
mtext(side=2,line=1,"Proportion of SNPs",outer=T,cex=1.2)

dev.off()

