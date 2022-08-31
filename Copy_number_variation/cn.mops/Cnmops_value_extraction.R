require(openxlsx)
#MK_locass<-read.xlsx("LocalAssessments_data_cnmops_MKar.xlsx",1)
MK_indCall<-read.xlsx("IndividualCall_data_cnmops_MKar.xlsx",1)
#CNV call produced by the segmentation algorithm
MK_normcov<-read.xlsx("Normalized_data_cnmops_MKar.xlsx",1)
#MK_posteriorprobs<-read.xlsx("Posteriorprobs_data_cnmops_MKar.xlsx",1)

#KM_locass<-read.xlsx("LocalAssessments_data_cnmops_KMar.xlsx",1)
KM_indCall<-read.xlsx("IndividualCall_data_cnmops_KMar.xlsx",1)
KM_normcov<-read.xlsx("Normalized_data_cnmops_KMar.xlsx",1)
#KM_posteriorprobs<-read.xlsx("Posteriorprobs_data_cnmops_KMar.xlsx",1)

#mean_KM_locass<-apply(KM_locass,1,mean)
mean_KM_indCall<-apply(KM_indCall,1,mean)
mean_KM_normcov<-apply(KM_normcov,1,mean)
#for (i in 1:8)
#	{assign(paste0("mean_KM_posteriorprobs_CN",i),apply(-log(KM_posteriorprobs[,c(1:14)+(14*(i-1))]+1),1,mean))
#	}

#mean_MK_locass<-apply(MK_locass,1,mean)
mean_MK_indCall<-apply(MK_indCall,1,mean)
mean_MK_normcov<-apply(MK_normcov,1,mean)
#for (i in 1:8)
#	{assign(paste0("mean_MK_posteriorprobs_CN",i),apply(-log(MK_posteriorprobs[,c(1:14)+(14*(i-1))]+1),1,mean))
#	}


MK<-data.frame(mean_MK_indCall,mean_MK_normcov)
KM<-data.frame(mean_KM_indCall,mean_KM_normcov)

Pos<-read.xlsx("Pos__data_cnmops_MKar.xlsx",1)
require(tidyr)
Pos<-separate(Pos,1,sep="_",into=c("NA","Scaff","Start","End"),remove=T)
Pos$Scaffold<-paste(Pos[,1],Pos[,2],sep="_")
Pos<-Pos[,c(5,3,4)]
Pos[,2]<-as.numeric(Pos[,2])
Pos[,3]<-as.numeric(Pos[,3])

MK2<-data.frame(Pos,MK)
KM2<-data.frame(Pos,KM)

require(dichromat)
colfunc <- colorRampPalette(c("black","red"))
colours_10<-colfunc(10)


pdf("Cnmops_values.pdf",width=8, height=12,paper="special")
par(mar=c(6,6,0,1)+0.1)
par(mfrow=c(3,1))
plot(mean_KM_locass,mean_KM_indCall)
plot(mean_KM_normcov,mean_KM_indCall)
plot(mean_KM_posteriorprobs_CN1,mean_KM_indCall,col=colours_10[1])
for (i in 2:8)
	{points(get(paste0("mean_KM_posteriorprobs_CN",i)),mean_KM_indCall,col=colours_10[i])}
dev.off()

pdf("Cnmops_values_hist.pdf",width=12, height=12,paper="special")
par(mar=c(6,6,2,1)+0.1)
par(mfrow=c(2,2))
hist(mean_KM_indCall)
hist(mean_KM_locass)
hist(mean_KM_normcov)
hist(mean_KM_posteriorprobs_CN1,col=colours_10[1])
for (i in 2:8)
	{hist(get(paste0("mean_KM_posteriorprobs_CN",i)),col=colours_10[i],add=T)}
dev.off()
options(scipen=999)
require(Hmisc)
pdf("Hist_Cnmops_values_MiasKowa_paper_square.pdf",width=10,height=8,paper="special",pointsize=25)
par(oma=c(1,2,1,0))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
hist(mean_MK_indCall,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,200000),xlim=c(-1.85,2))
mgp.axis(pos=-2,side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=axTicks(side=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"CNV call",line=-2,outer=T,cex=1.5)

hist(mean_MK_normcov,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,600000))
mgp.axis(pos=min(hist(mean_MK_normcov,breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=formatC(axTicks(side=1),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"Coverage",line=-2,outer=T,cex=1.5)

dev.off()


pdf("Hist_Cnmops_values_MiasKowa_paper_square_log2.pdf",width=10,height=7.5,paper="special",pointsize=28)
par(oma=c(1,2,1,0))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
options(scipen=999)
hist(log10(mean_MK_normcov),breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,100000))
mgp.axis(pos=min(hist(log10(mean_MK_normcov),breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2)/1000,big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1)[c(1,3,5,7)],labels=c(formatC(10^(axTicks(side=1)[c(1,3,5)]),big.mark=",",format="g",decimal.mark = "."),formatC(10^(axTicks(side=1)[c(7)]),big.mark=",",format="g",decimal.mark = ".",digits=6)),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=0,outer=T,cex=1.5)
mtext(side=1,"Coverage",line=-2,outer=T,cex=1.5)

dev.off()


pdf("Hist_Cnmops_values_MiasKowa_paper_square_log.pdf",width=10,height=8,paper="special",pointsize=28)
par(oma=c(1,2,1,0))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
hist(mean_MK_indCall,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,200000),xlim=c(-1.85,2))
mgp.axis(pos=-2,side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2)/1000,big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=axTicks(side=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"CNV call",line=-2,outer=T,cex=1.5)

options(scipen=999)
hist(log10(mean_MK_normcov),breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,100000))
mgp.axis(pos=min(hist(log10(mean_MK_normcov),breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2)/1000,big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1)[c(1,3,5,7)],labels=c(formatC(10^(axTicks(side=1)[c(1,3,5)]),big.mark=",",format="g",decimal.mark = "."),formatC(10^(axTicks(side=1)[c(7)]),big.mark=",",format="g",decimal.mark = ".",digits=6)),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=0,outer=T,cex=1.5)
mtext(side=1,"Coverage",line=-2,outer=T,cex=1.5)

dev.off()

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


#HMA4 AL3G52820 scaffold_3 23475941 23483785 -
scaffold="scaffold_3"
start=23475941
stop=23483785
windowsize=50000
HMA4_MK<-MK2[MK2$Scaffold==scaffold&MK2$Start>=(start-windowsize)&MK2$End<=(stop+windowsize),]
HMA4_KM<-KM2[KM2$Scaffold==scaffold&KM2$Start>=(start-windowsize)&KM2$End<=(stop+windowsize),]
HMA4_MK$Mid<-HMA4_MK$Start+((HMA4_MK$End-HMA4_MK$Start)/2)
HMA4_KM$Mid<-HMA4_KM$Start+((HMA4_MK$End-HMA4_MK$Start)/2)

pdf("Cnmops_Mias_HMA4_paper_test.pdf",width=10,height=15,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G52820",]
par(mar=c(1,6,0,1))
par(mgp=c(2.5,0.75,0))
par(oma=c(5,2,2,2))

arrowdir2<-1#-
par(mfrow=c(3,1))
#mgp.axis(side=2,at=axTicks(side=2),labels=c(axTicks(side=2)[-length(axTicks(side=2))],NA),cex=1,cex.lab=1.2,cex.axis=1.2,las=1)
plot(HMA4_MK$mean_MK_locass~HMA4_MK$Mid,col="red",type="l",ylab="Local Assessments",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,xlab="",ylim=c(min(HMA4_MK$mean_MK_locass,HMA4_KM$mean_KM_locass),max(HMA4_MK$mean_MK_locass+2,HMA4_KM$mean_KM_locass+2)))
rect(start,-100,stop,100,col="lightgrey",border = NA)	
lines(HMA4_MK$mean_MK_locass~HMA4_MK$Mid,col="red",lwd=2)
lines(HMA4_KM$mean_KM_locass~HMA4_KM$Mid,col="black",lwd=2)
box()

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=Lyr_Genesofinterest$gene_end[i],y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=arrowdir[i],length=0.1,col="grey40",lwd=1)
	}
arrows(x0=start,y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=stop,y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=Exonsofinterest$gene_end[i],y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=Exonsofinterestcand$gene_end[i],y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((Start-End)/2),1.2,"HMA4",col="red",font=3)
#legend("topleft",lwd=2,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)

plot(HMA4_MK$mean_MK_indCall~HMA4_MK$Mid,col="red",type="l",xlab="",ylab="CN call",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK$mean_MK_indCall,HMA4_KM$mean_KM_indCall),max(HMA4_MK$mean_MK_indCall,HMA4_KM$mean_KM_indCall)))
rect(start,-100,stop,100,col="lightgrey",border = NA)	
lines(HMA4_MK$mean_MK_indCall~HMA4_MK$Mid,col="red",lwd=2)
lines(HMA4_KM$mean_KM_indCall~HMA4_KM$Mid,col="black",lwd=2)
box()

plot(HMA4_MK$mean_MK_normcov~HMA4_MK$Mid,col="red",type="l",xlab="Scaffold position (Mb)",ylab="Normalized coverage",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK$mean_MK_normcov,HMA4_KM$mean_KM_normcov),max(HMA4_MK$mean_MK_normcov,HMA4_KM$mean_KM_normcov)))
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="g",decimal.mark = "."),cex=1,cex.lab=1.2,cex.axis=1.2,xpd=T)
lines(HMA4_MK$mean_MK_normcov~HMA4_MK$Mid,col="red",lwd=2)
lines(HMA4_KM$mean_KM_normcov~HMA4_KM$Mid,col="black",lwd=2)
box()
dev.off()

require(Hmisc)
#HMA4
#scaffold_3	23475941 23483785	9845
pdf("Cnmops_Mias_HMA4_paper.pdf",width=10,height=15,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G52820",]
par(mar=c(1,6,0,1))
par(mgp=c(2.5,0.75,0))
par(oma=c(5,2,2,2))

arrowdir2<-1#-
par(mfrow=c(3,1))
#mgp.axis(side=2,at=axTicks(side=2),labels=c(axTicks(side=2)[-length(axTicks(side=2))],NA),cex=1,cex.lab=1.2,cex.axis=1.2,las=1)
plot(HMA4_MK$mean_MK_locass~HMA4_MK$Mid,col="red",type="l",ylab="Local Assessments",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,xlab="",ylim=c(min(HMA4_MK$mean_MK_locass,HMA4_KM$mean_KM_locass),max(HMA4_MK$mean_MK_locass+2,HMA4_KM$mean_KM_locass+2)))
rect(start,-100,stop,100,col="lightgrey",border = NA)	
lines(HMA4_MK$mean_MK_locass~HMA4_MK$Mid,col="red",lwd=2)
lines(HMA4_KM$mean_KM_locass~HMA4_KM$Mid,col="black",lwd=2)
box()

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=Lyr_Genesofinterest$gene_end[i],y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=arrowdir[i],length=0.1,col="grey40",lwd=1)
	}
arrows(x0=start,y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=stop,y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=Exonsofinterest$gene_end[i],y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),x1=Exonsofinterestcand$gene_end[i],y1=max(HMA4_MK$mean_MK_locass+1,HMA4_KM$mean_KM_locass+1),code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((Start-End)/2),1.2,"HMA4",col="red",font=3)
#legend("topleft",lwd=2,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)

plot(HMA4_MK$mean_MK_indCall~HMA4_MK$Mid,col="red",type="l",xlab="",ylab="CN call",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK$mean_MK_indCall,HMA4_KM$mean_KM_indCall),max(HMA4_MK$mean_MK_indCall,HMA4_KM$mean_KM_indCall)))
rect(start,-100,stop,100,col="lightgrey",border = NA)	
lines(HMA4_MK$mean_MK_indCall~HMA4_MK$Mid,col="red",lwd=2)
lines(HMA4_KM$mean_KM_indCall~HMA4_KM$Mid,col="black",lwd=2)
box()

plot(HMA4_MK$mean_MK_normcov~HMA4_MK$Mid,col="red",type="l",xlab="Scaffold position (Mb)",ylab="Normalized coverage",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK$mean_MK_normcov,HMA4_KM$mean_KM_normcov),max(HMA4_MK$mean_MK_normcov,HMA4_KM$mean_KM_normcov)))
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="g",decimal.mark = "."),cex=1,cex.lab=1.2,cex.axis=1.2,xpd=T)
lines(HMA4_MK$mean_MK_normcov~HMA4_MK$Mid,col="red",lwd=2)
lines(HMA4_KM$mean_KM_normcov~HMA4_KM$Mid,col="black",lwd=2)
box()
dev.off()


#MK_locass2<-data.frame(Pos,MK_locass)
#KM_locass2<-data.frame(Pos,KM_locass)
MK_indCall2<-data.frame(Pos,MK_indCall)
KM_indCall2<-data.frame(Pos,KM_indCall)
MK_normcov2<-data.frame(Pos,MK_normcov)
KM_normcov2<-data.frame(Pos,KM_normcov)

#HMA4_MK_locass<-MK_locass2[MK_locass2$Scaffold==scaffold&MK_locass2$Start>=(start-windowsize)&MK_locass2$End<=(stop+windowsize),]
#HMA4_KM_locass<-KM_locass2[KM_locass2$Scaffold==scaffold&KM_locass2$Start>=(start-windowsize)&KM_locass2$End<=(stop+windowsize),]
#HMA4_MK_locass$Mid<-HMA4_MK_locass$Start+((HMA4_MK_locass$End-HMA4_MK_locass$Start)/2)
#HMA4_KM_locass$Mid<-HMA4_KM_locass$Start+((HMA4_MK_locass$End-HMA4_MK_locass$Start)/2)

HMA4_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
HMA4_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
HMA4_MK_indCall$Mid<-HMA4_MK_indCall$Start+((HMA4_MK_indCall$End-HMA4_MK_indCall$Start)/2)
HMA4_KM_indCall$Mid<-HMA4_KM_indCall$Start+((HMA4_MK_indCall$End-HMA4_MK_indCall$Start)/2)

HMA4_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
HMA4_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
HMA4_MK_normcov$Mid<-HMA4_MK_normcov$Start+((HMA4_MK_normcov$End-HMA4_MK_normcov$Start)/2)
HMA4_KM_normcov$Mid<-HMA4_KM_normcov$Start+((HMA4_MK_normcov$End-HMA4_MK_normcov$Start)/2)


pdf("Cnmops_Mias_HMA4_paper_inds_test.pdf",width=10,height=15,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G52820",]
par(mar=c(1,6,0,1))
par(mgp=c(2.5,0.75,0))
par(oma=c(5,2,2,2))

arrowdir2<-1#-
par(mfrow=c(3,1))
#mgp.axis(side=2,at=axTicks(side=2),labels=c(axTicks(side=2)[-length(axTicks(side=2))],NA),cex=1,cex.lab=1.2,cex.axis=1.2,las=1)
plot(HMA4_MK_locass[,4]~HMA4_MK_locass$Mid,col="red",type="l",ylab="Local Assessments",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,xlab="",ylim=c(min(HMA4_MK_locass[,4:22],HMA4_KM_locass[,4:11]),max(HMA4_MK_locass[,4:22]+2,HMA4_KM_locass[,4:11]+2)))
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(HMA4_MK_locass[,i]~HMA4_MK_locass$Mid,col="red",lwd=2)}
for (j in 4:11)
	{lines(HMA4_KM_locass[,j]~HMA4_KM_locass$Mid,col="black",lwd=2)}
box()

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),x1=Lyr_Genesofinterest$gene_end[i],y1=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),code=arrowdir[i],length=0.1,col="grey40",lwd=1)
	}
arrows(x0=start,y0=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),x1=stop,y1=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),x1=Exonsofinterest$gene_end[i],y1=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),x1=Exonsofinterestcand$gene_end[i],y1=max(HMA4_MK_locass[,4:22]+1,HMA4_KM_locass[,4:11]+1),code=arrowdir2,length=0,col="red",lwd=5)
	}
#text(Start+((Start-End)/2),1.2,"HMA4",col="red",font=3)
#legend("topleft",lwd=2,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T, inset=c(-0.1,-0.2),xpd=T)

plot(HMA4_MK_indCall[,4]~HMA4_MK_indCall$Mid,col="red",type="l",xlab="",ylab="CN call",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK_indCall[,4:22],HMA4_KM_indCall[,4:11]),max(HMA4_MK_indCall[,4:22],HMA4_KM_indCall[,4:11])))
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(HMA4_MK_indCall[,i]~HMA4_MK_indCall$Mid,col="red",lwd=2)}
for (j in 4:11)
	{lines(HMA4_KM_indCall[,j]~HMA4_KM_indCall$Mid,col="black",lwd=2)}
box()

plot(HMA4_MK_normcov[,4]~HMA4_MK_normcov$Mid,col="red",type="l",xlab="Scaffold position (Mb)",ylab="Normalized coverage",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK_normcov[,4:22],HMA4_KM_normcov[,4:11]),max(HMA4_MK_normcov[,4:22],HMA4_KM_normcov[,4:11])))
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="g",decimal.mark = "."),cex=1,cex.lab=1.2,cex.axis=1.2,xpd=T)
for (i in 4:22)
	{lines(HMA4_MK_normcov[,i]~HMA4_MK_normcov$Mid,col="red",lwd=2)}
for (j in 4:11)
	{lines(HMA4_KM_normcov[,j]~HMA4_KM_normcov$Mid,col="black",lwd=2)}
box()
dev.off()


require(Hmisc)
#HMA4
#scaffold_3	23475941 23483785	9845
pdf("Cnmops_Mias_HMA4_paper_inds.pdf",width=15,height=10,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G52820",]
par(mar=c(0,6,2,1))
par(mgp=c(3,0.75,0))
par(oma=c(5,2,2,2))

arrowdir2<-1#-
par(mfrow=c(2,1))

plot(HMA4_MK_indCall[,4]~HMA4_MK_indCall$Mid,col="red",type="l",xlab="",ylab="CN I/NI call",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK_indCall[,4:22],HMA4_KM_indCall[,4:11]),max(HMA4_MK_indCall[,4:22]+2,HMA4_KM_indCall[,4:11]+2)))
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(HMA4_MK_indCall[,i]~HMA4_MK_indCall$Mid,col="red",lwd=2)}
for (j in 4:11)
	{lines(HMA4_KM_indCall[,j]~HMA4_KM_indCall$Mid,col="black",lwd=2)}
box()
legend("topleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.2),xpd=TRUE)

arrowdir<-ifelse(Lyr_Genesofinterest$Strand=="+",2,1)
for(i in 1:nrow(Lyr_Genesofinterest)) 
	{arrows(x0=Lyr_Genesofinterest$gene_start[i],y0=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),x1=Lyr_Genesofinterest$gene_end[i],y1=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),code=arrowdir[i],length=0.1,col="grey40",lwd=1)
	}
arrows(x0=start,y0=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),x1=stop,y1=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),code=arrowdir2,length=0.1,col="red",lwd=1)
for(i in 1:nrow(Exonsofinterest)) 
	{arrows(x0=Exonsofinterest$gene_start[i],y0=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),x1=Exonsofinterest$gene_end[i],y1=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),code=1,length=0,col="grey40",lwd=5)
	}
for(i in 1:nrow(Exonsofinterestcand)) 
	{arrows(x0=Exonsofinterestcand$gene_start[i],y0=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),x1=Exonsofinterestcand$gene_end[i],y1=max(HMA4_MK_indCall[,4:22]+1,HMA4_KM_indCall[,4:11]+1),code=arrowdir2,length=0,col="red",lwd=5)
	}
plot(HMA4_MK_normcov[,4]~HMA4_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK_normcov[,4:22],HMA4_KM_normcov[,4:11]),max(HMA4_MK_normcov[,4:22],HMA4_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="g",decimal.mark = "."),cex=1,cex.lab=1.2,cex.axis=1.2,xpd=T)
for (i in 4:22)
	{lines(HMA4_MK_normcov[,i]~HMA4_MK_normcov$Mid,col="red",lwd=2)}
for (j in 4:11)
	{lines(HMA4_KM_normcov[,j]~HMA4_KM_normcov$Mid,col="black",lwd=2)}
box()
mtext(side=2,"Normalized coverage",line=3,xpd=T,cex=1.5,adj=0.3)
mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
dev.off()

pdf("Cnmops_Mias_HMA4_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G52820",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-1#-
par(mfrow=c(2,1))

plot(HMA4_MK_indCall[,4]~HMA4_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK_indCall[,4:22],HMA4_KM_indCall[,4:11]),max(HMA4_MK_indCall[,4:22],HMA4_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(HMA4_MK_indCall[,i]~HMA4_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA4_KM_indCall[,j]~HMA4_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA4_MK_normcov[,4]~HMA4_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_MK_normcov[,4:22],HMA4_KM_normcov[,4:11]),max(HMA4_MK_normcov[,4:22],HMA4_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(HMA4_MK_normcov[,i]~HMA4_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA4_KM_normcov[,j]~HMA4_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#MTP1
scaffold="scaffold_4"
start=22776960
stop=22778882
windowsize=50000

MTP1_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
MTP1_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
MTP1_MK_indCall$Mid<-MTP1_MK_indCall$Start+((MTP1_MK_indCall$End-MTP1_MK_indCall$Start)/2)
MTP1_KM_indCall$Mid<-MTP1_KM_indCall$Start+((MTP1_MK_indCall$End-MTP1_MK_indCall$Start)/2)

MTP1_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
MTP1_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
MTP1_MK_normcov$Mid<-MTP1_MK_normcov$Start+((MTP1_MK_normcov$End-MTP1_MK_normcov$Start)/2)
MTP1_KM_normcov$Mid<-MTP1_KM_normcov$Start+((MTP1_MK_normcov$End-MTP1_MK_normcov$Start)/2)

pdf("Cnmops_Mias_MTP1_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL4G46270",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-2#+
par(mfrow=c(2,1))

plot(MTP1_MK_indCall[,4]~MTP1_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(MTP1_MK_indCall[,4:22],MTP1_KM_indCall[,4:11]),max(MTP1_MK_indCall[,4:22],MTP1_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(MTP1_MK_indCall[,i]~MTP1_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(MTP1_KM_indCall[,j]~MTP1_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(MTP1_MK_normcov[,4]~MTP1_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(MTP1_MK_normcov[,4:22],MTP1_KM_normcov[,4:11]),max(MTP1_MK_normcov[,4:22],MTP1_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(MTP1_MK_normcov[,i]~MTP1_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(MTP1_KM_normcov[,j]~MTP1_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


#HMA3
scaffold="scaffold_7"
start=4960535
stop=4963880
windowsize=50000

HMA3_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
HMA3_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
HMA3_MK_indCall$Mid<-HMA3_MK_indCall$Start+((HMA3_MK_indCall$End-HMA3_MK_indCall$Start)/2)
HMA3_KM_indCall$Mid<-HMA3_KM_indCall$Start+((HMA3_MK_indCall$End-HMA3_MK_indCall$Start)/2)

HMA3_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
HMA3_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
HMA3_MK_normcov$Mid<-HMA3_MK_normcov$Start+((HMA3_MK_normcov$End-HMA3_MK_normcov$Start)/2)
HMA3_KM_normcov$Mid<-HMA3_KM_normcov$Start+((HMA3_MK_normcov$End-HMA3_MK_normcov$Start)/2)

pdf("Cnmops_Mias_HMA3_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G22080",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-2#+
par(mfrow=c(2,1))

plot(HMA3_MK_indCall[,4]~HMA3_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA3_MK_indCall[,4:22],HMA3_KM_indCall[,4:11]),max(HMA3_MK_indCall[,4:22],HMA3_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(HMA3_MK_indCall[,i]~HMA3_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA3_KM_indCall[,j]~HMA3_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA3_MK_normcov[,4]~HMA3_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA3_MK_normcov[,4:22],HMA3_KM_normcov[,4:11]),max(HMA3_MK_normcov[,4:22],HMA3_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(HMA3_MK_normcov[,i]~HMA3_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA3_KM_normcov[,j]~HMA3_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP6
scaffold="scaffold_4"
start=13949254
stop=13950804
windowsize=50000

ZIP6_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
ZIP6_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
ZIP6_MK_indCall$Mid<-ZIP6_MK_indCall$Start+((ZIP6_MK_indCall$End-ZIP6_MK_indCall$Start)/2)
ZIP6_KM_indCall$Mid<-ZIP6_KM_indCall$Start+((ZIP6_MK_indCall$End-ZIP6_MK_indCall$Start)/2)

ZIP6_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
ZIP6_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
ZIP6_MK_normcov$Mid<-ZIP6_MK_normcov$Start+((ZIP6_MK_normcov$End-ZIP6_MK_normcov$Start)/2)
ZIP6_KM_normcov$Mid<-ZIP6_KM_normcov$Start+((ZIP6_MK_normcov$End-ZIP6_MK_normcov$Start)/2)

pdf("Cnmops_Mias_ZIP6_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL4G46270",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-1#-
par(mfrow=c(2,1))

plot(ZIP6_MK_indCall[,4]~ZIP6_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP6_MK_indCall[,4:22],ZIP6_KM_indCall[,4:11]),max(ZIP6_MK_indCall[,4:22],ZIP6_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(ZIP6_MK_indCall[,i]~ZIP6_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP6_KM_indCall[,j]~ZIP6_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP6_MK_normcov[,4]~ZIP6_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP6_MK_normcov[,4:22],ZIP6_KM_normcov[,4:11]),max(ZIP6_MK_normcov[,4:22],ZIP6_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(ZIP6_MK_normcov[,i]~ZIP6_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP6_KM_normcov[,j]~ZIP6_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#NRAMP3
scaffold="scaffold_4"
start=1634727
stop=1637073
windowsize=50000

NRAMP3_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
NRAMP3_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
NRAMP3_MK_indCall$Mid<-NRAMP3_MK_indCall$Start+((NRAMP3_MK_indCall$End-NRAMP3_MK_indCall$Start)/2)
NRAMP3_KM_indCall$Mid<-NRAMP3_KM_indCall$Start+((NRAMP3_MK_indCall$End-NRAMP3_MK_indCall$Start)/2)

NRAMP3_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
NRAMP3_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
NRAMP3_MK_normcov$Mid<-NRAMP3_MK_normcov$Start+((NRAMP3_MK_normcov$End-NRAMP3_MK_normcov$Start)/2)
NRAMP3_KM_normcov$Mid<-NRAMP3_KM_normcov$Start+((NRAMP3_MK_normcov$End-NRAMP3_MK_normcov$Start)/2)

pdf("Cnmops_Mias_NRAMP3_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL4G13010",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-1#-
par(mfrow=c(2,1))

plot(NRAMP3_MK_indCall[,4]~NRAMP3_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(NRAMP3_MK_indCall[,4:22],NRAMP3_KM_indCall[,4:11]),max(NRAMP3_MK_indCall[,4:22],NRAMP3_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(NRAMP3_MK_indCall[,i]~NRAMP3_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(NRAMP3_KM_indCall[,j]~NRAMP3_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(NRAMP3_MK_normcov[,4]~NRAMP3_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(NRAMP3_MK_normcov[,4:22],NRAMP3_KM_normcov[,4:11]),max(NRAMP3_MK_normcov[,4:22],NRAMP3_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,2000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(NRAMP3_MK_normcov[,i]~NRAMP3_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(NRAMP3_KM_normcov[,j]~NRAMP3_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP11
scaffold="scaffold_1"
start=29573151
stop=29574444
windowsize=50000

ZIP11_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
ZIP11_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
ZIP11_MK_indCall$Mid<-ZIP11_MK_indCall$Start+((ZIP11_MK_indCall$End-ZIP11_MK_indCall$Start)/2)
ZIP11_KM_indCall$Mid<-ZIP11_KM_indCall$Start+((ZIP11_MK_indCall$End-ZIP11_MK_indCall$Start)/2)

ZIP11_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
ZIP11_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
ZIP11_MK_normcov$Mid<-ZIP11_MK_normcov$Start+((ZIP11_MK_normcov$End-ZIP11_MK_normcov$Start)/2)
ZIP11_KM_normcov$Mid<-ZIP11_KM_normcov$Start+((ZIP11_MK_normcov$End-ZIP11_MK_normcov$Start)/2)

pdf("Cnmops_Mias_ZIP11_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL1G64050",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-2#+
par(mfrow=c(2,1))

plot(ZIP11_MK_indCall[,4]~ZIP11_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP11_MK_indCall[,4:22],ZIP11_KM_indCall[,4:11]),max(ZIP11_MK_indCall[,4:22],ZIP11_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(ZIP11_MK_indCall[,i]~ZIP11_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP11_KM_indCall[,j]~ZIP11_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP11_MK_normcov[,4]~ZIP11_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP11_MK_normcov[,4:22],ZIP11_KM_normcov[,4:11]),max(ZIP11_MK_normcov[,4:22],ZIP11_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(ZIP11_MK_normcov[,i]~ZIP11_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP11_KM_normcov[,j]~ZIP11_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#HMA2
scaffold="scaffold_7"
start=4971705
stop=4976723
windowsize=50000

HMA2_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
HMA2_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
HMA2_MK_indCall$Mid<-HMA2_MK_indCall$Start+((HMA2_MK_indCall$End-HMA2_MK_indCall$Start)/2)
HMA2_KM_indCall$Mid<-HMA2_KM_indCall$Start+((HMA2_MK_indCall$End-HMA2_MK_indCall$Start)/2)

HMA2_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
HMA2_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
HMA2_MK_normcov$Mid<-HMA2_MK_normcov$Start+((HMA2_MK_normcov$End-HMA2_MK_normcov$Start)/2)
HMA2_KM_normcov$Mid<-HMA2_KM_normcov$Start+((HMA2_MK_normcov$End-HMA2_MK_normcov$Start)/2)

pdf("Cnmops_Mias_HMA2_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G22090",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-2#+
par(mfrow=c(2,1))

plot(HMA2_MK_indCall[,4]~HMA2_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA2_MK_indCall[,4:22],HMA2_KM_indCall[,4:11]),max(HMA2_MK_indCall[,4:22],HMA2_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(HMA2_MK_indCall[,i]~HMA2_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA2_KM_indCall[,j]~HMA2_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA2_MK_normcov[,4]~HMA2_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA2_MK_normcov[,4:22],HMA2_KM_normcov[,4:11]),max(HMA2_MK_normcov[,4:22],HMA2_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(HMA2_MK_normcov[,i]~HMA2_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA2_KM_normcov[,j]~HMA2_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#IRT2
scaffold="scaffold_7"
start=10093506
stop=10094856
windowsize=50000

IRT2_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
IRT2_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
IRT2_MK_indCall$Mid<-IRT2_MK_indCall$Start+((IRT2_MK_indCall$End-IRT2_MK_indCall$Start)/2)
IRT2_KM_indCall$Mid<-IRT2_KM_indCall$Start+((IRT2_MK_indCall$End-IRT2_MK_indCall$Start)/2)

IRT2_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
IRT2_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
IRT2_MK_normcov$Mid<-IRT2_MK_normcov$Start+((IRT2_MK_normcov$End-IRT2_MK_normcov$Start)/2)
IRT2_KM_normcov$Mid<-IRT2_KM_normcov$Start+((IRT2_MK_normcov$End-IRT2_MK_normcov$Start)/2)

pdf("Cnmops_Mias_IRT2_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL7G34360",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-2#+
par(mfrow=c(2,1))

plot(IRT2_MK_indCall[,4]~IRT2_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(IRT2_MK_indCall[,4:22],IRT2_KM_indCall[,4:11]),max(IRT2_MK_indCall[,4:22],IRT2_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(IRT2_MK_indCall[,i]~IRT2_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(IRT2_KM_indCall[,j]~IRT2_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(IRT2_MK_normcov[,4]~IRT2_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(IRT2_MK_normcov[,4:22],IRT2_KM_normcov[,4:11]),max(IRT2_MK_normcov[,4:22],IRT2_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(IRT2_MK_normcov[,i]~IRT2_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(IRT2_KM_normcov[,j]~IRT2_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP8
scaffold="scaffold_8"
start=2517313
stop=2519529
windowsize=50000

ZIP8_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
ZIP8_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
ZIP8_MK_indCall$Mid<-ZIP8_MK_indCall$Start+((ZIP8_MK_indCall$End-ZIP8_MK_indCall$Start)/2)
ZIP8_KM_indCall$Mid<-ZIP8_KM_indCall$Start+((ZIP8_MK_indCall$End-ZIP8_MK_indCall$Start)/2)

ZIP8_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
ZIP8_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
ZIP8_MK_normcov$Mid<-ZIP8_MK_normcov$Start+((ZIP8_MK_normcov$End-ZIP8_MK_normcov$Start)/2)
ZIP8_KM_normcov$Mid<-ZIP8_KM_normcov$Start+((ZIP8_MK_normcov$End-ZIP8_MK_normcov$Start)/2)

pdf("Cnmops_Mias_ZIP8_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL8G14120",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-2#+
par(mfrow=c(2,1))

plot(ZIP8_MK_indCall[,4]~ZIP8_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP8_MK_indCall[,4:22],ZIP8_KM_indCall[,4:11]),max(ZIP8_MK_indCall[,4:22],ZIP8_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(ZIP8_MK_indCall[,i]~ZIP8_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP8_KM_indCall[,j]~ZIP8_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP8_MK_normcov[,4]~ZIP8_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP8_MK_normcov[,4:22],ZIP8_KM_normcov[,4:11]),max(ZIP8_MK_normcov[,4:22],ZIP8_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(ZIP8_MK_normcov[,i]~ZIP8_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP8_KM_normcov[,j]~ZIP8_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


#FRO5
scaffold="scaffold_6"
start=10296305
stop=10299581
windowsize=50000

FRO5_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
FRO5_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
FRO5_MK_indCall$Mid<-FRO5_MK_indCall$Start+((FRO5_MK_indCall$End-FRO5_MK_indCall$Start)/2)
FRO5_KM_indCall$Mid<-FRO5_KM_indCall$Start+((FRO5_MK_indCall$End-FRO5_MK_indCall$Start)/2)

FRO5_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
FRO5_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
FRO5_MK_normcov$Mid<-FRO5_MK_normcov$Start+((FRO5_MK_normcov$End-FRO5_MK_normcov$Start)/2)
FRO5_KM_normcov$Mid<-FRO5_KM_normcov$Start+((FRO5_MK_normcov$End-FRO5_MK_normcov$Start)/2)

pdf("Cnmops_Mias_FRO5_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL6G35310",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-1#-
par(mfrow=c(2,1))

plot(FRO5_MK_indCall[,4]~FRO5_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FRO5_MK_indCall[,4:22],FRO5_KM_indCall[,4:11]),max(FRO5_MK_indCall[,4:22],FRO5_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(FRO5_MK_indCall[,i]~FRO5_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRO5_KM_indCall[,j]~FRO5_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FRO5_MK_normcov[,4]~FRO5_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FRO5_MK_normcov[,4:22],FRO5_KM_normcov[,4:11]),max(FRO5_MK_normcov[,4:22],FRO5_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(FRO5_MK_normcov[,i]~FRO5_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRO5_KM_normcov[,j]~FRO5_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#FER1
scaffold="scaffold_6"
start=187631
stop=189383
windowsize=50000

FER1_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
FER1_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
FER1_MK_indCall$Mid<-FER1_MK_indCall$Start+((FER1_MK_indCall$End-FER1_MK_indCall$Start)/2)
FER1_KM_indCall$Mid<-FER1_KM_indCall$Start+((FER1_MK_indCall$End-FER1_MK_indCall$Start)/2)

FER1_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
FER1_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
FER1_MK_normcov$Mid<-FER1_MK_normcov$Start+((FER1_MK_normcov$End-FER1_MK_normcov$Start)/2)
FER1_KM_normcov$Mid<-FER1_KM_normcov$Start+((FER1_MK_normcov$End-FER1_MK_normcov$Start)/2)

pdf("Cnmops_Mias_FER1_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL6G10450",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-2#+
par(mfrow=c(2,1))

plot(FER1_MK_indCall[,4]~FER1_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FER1_MK_indCall[,4:22],FER1_KM_indCall[,4:11]),max(FER1_MK_indCall[,4:22],FER1_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(FER1_MK_indCall[,i]~FER1_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FER1_KM_indCall[,j]~FER1_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FER1_MK_normcov[,4]~FER1_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FER1_MK_normcov[,4:22],FER1_KM_normcov[,4:11]),max(FER1_MK_normcov[,4:22],FER1_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(FER1_MK_normcov[,i]~FER1_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FER1_KM_normcov[,j]~FER1_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#FRD3
scaffold="scaffold_3"
start=3352004
stop=3357699
windowsize=50000

FRD3_MK_indCall<-MK_indCall2[MK_indCall2$Scaffold==scaffold&MK_indCall2$Start>=(start-windowsize)&MK_indCall2$End<=(stop+windowsize),]
FRD3_KM_indCall<-KM_indCall2[KM_indCall2$Scaffold==scaffold&KM_indCall2$Start>=(start-windowsize)&KM_indCall2$End<=(stop+windowsize),]
FRD3_MK_indCall$Mid<-FRD3_MK_indCall$Start+((FRD3_MK_indCall$End-FRD3_MK_indCall$Start)/2)
FRD3_KM_indCall$Mid<-FRD3_KM_indCall$Start+((FRD3_MK_indCall$End-FRD3_MK_indCall$Start)/2)

FRD3_MK_normcov<-MK_normcov2[MK_normcov2$Scaffold==scaffold&MK_normcov2$Start>=(start-windowsize)&MK_normcov2$End<=(stop+windowsize),]
FRD3_KM_normcov<-KM_normcov2[KM_normcov2$Scaffold==scaffold&KM_normcov2$Start>=(start-windowsize)&KM_normcov2$End<=(stop+windowsize),]
FRD3_MK_normcov$Mid<-FRD3_MK_normcov$Start+((FRD3_MK_normcov$End-FRD3_MK_normcov$Start)/2)
FRD3_KM_normcov$Mid<-FRD3_KM_normcov$Start+((FRD3_MK_normcov$End-FRD3_MK_normcov$Start)/2)

pdf("Cnmops_Mias_FRD3_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
Lyr_Genes1<-genelist[genelist$Scaffold==scaffold,]
Lyr_Genesofinterest<-Lyr_Genes1[Lyr_Genes1$gene_start>=(start-windowsize)&Lyr_Genes1$gene_start<=(start+windowsize)|Lyr_Genes1$gene_end>=(stop-windowsize)&Lyr_Genes1$gene_end<=(stop+windowsize),]
Exonsofinterest1<-exonlist[exonlist$Scaffold==scaffold,]
Exonsofinterest<-Exonsofinterest1[Exonsofinterest1$gene_start>=(start-windowsize)&Exonsofinterest1$gene_start<=(start+windowsize)|Exonsofinterest1$gene_end>=(stop-windowsize)&Exonsofinterest1$gene_end<=(stop+windowsize),]
Exonsofinterestcand<-exonlist[exonlist$Lyr_Gene=="AL3G19430",]
par(mar=c(0,6,1,1)+0.1)
par(oma=c(1,1,1,0))
par(mgp=c(2.5,0.75,0))
par(lwd=2.5)

arrowdir2<-1#-
par(mfrow=c(2,1))

plot(FRD3_MK_indCall[,4]~FRD3_MK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FRD3_MK_indCall[,4:22],FRD3_KM_indCall[,4:11]),max(FRD3_MK_indCall[,4:22],FRD3_KM_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:22)
	{lines(FRD3_MK_indCall[,i]~FRD3_MK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRD3_KM_indCall[,j]~FRD3_KM_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FRD3_MK_normcov[,4]~FRD3_MK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FRD3_MK_normcov[,4:22],FRD3_KM_normcov[,4:11]),max(FRD3_MK_normcov[,4:22],FRD3_KM_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:22)
	{lines(FRD3_MK_normcov[,i]~FRD3_MK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRD3_KM_normcov[,j]~FRD3_KM_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Mias"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


