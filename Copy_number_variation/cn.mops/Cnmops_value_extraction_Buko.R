require(openxlsx)
BK_indCall<-read.xlsx("IndividualCall_data_cnmops_BKar.xlsx",1)
BK_normcov<-read.xlsx("Normalized_data_cnmops_BKar.xlsx",1)

KB_indCall<-read.xlsx("IndividualCall_data_cnmops_KBar.xlsx",1)
KB_normcov<-read.xlsx("Normalized_data_cnmops_KBar.xlsx",1)

mean_KB_indCall<-apply(KB_indCall,1,mean)
mean_KB_normcov<-apply(KB_normcov,1,mean)
mean_BK_indCall<-apply(BK_indCall,1,mean)
mean_BK_normcov<-apply(BK_normcov,1,mean)

BK<-data.frame(mean_BK_indCall,mean_BK_normcov)
KB<-data.frame(mean_KB_indCall,mean_KB_normcov)

Pos<-read.xlsx("Pos__data_cnmops_MKar.xlsx",1)
require(tidyr)
Pos<-separate(Pos,1,sep="_",into=c("NA","Scaff","Start","End"),remove=T)
Pos$Scaffold<-paste(Pos[,1],Pos[,2],sep="_")
Pos<-Pos[,c(5,3,4)]
Pos[,2]<-as.numeric(Pos[,2])
Pos[,3]<-as.numeric(Pos[,3])

BK2<-data.frame(Pos,BK)
KB2<-data.frame(Pos,KB)

require(dichromat)
colfunc <- colorRampPalette(c("black","red"))
colours_10<-colfunc(10)


pdf("Cnmops_values.pdf",width=8, height=12,paper="special")
par(mar=c(6,6,0,1)+0.1)
par(mfrow=c(3,1))
plot(mean_KB_locass,mean_KB_indCall)
plot(mean_KB_normcov,mean_KB_indCall)
plot(mean_KB_posteriorprobs_CN1,mean_KB_indCall,col=colours_10[1])
for (i in 2:8)
	{points(get(paste0("mean_KB_posteriorprobs_CN",i)),mean_KB_indCall,col=colours_10[i])}
dev.off()

pdf("Cnmops_values_hist.pdf",width=12, height=12,paper="special")
par(mar=c(6,6,2,1)+0.1)
par(mfrow=c(2,2))
hist(mean_KB_indCall)
hist(mean_KB_locass)
hist(mean_KB_normcov)
hist(mean_KB_posteriorprobs_CN1,col=colours_10[1])
for (i in 2:8)
	{hist(get(paste0("mean_KB_posteriorprobs_CN",i)),col=colours_10[i],add=T)}
dev.off()
options(scipen=999)
require(Hmisc)
pdf("Hist_Cnmops_values_BukoKowa_paper_square.pdf",width=10,height=8,paper="special",pointsize=25)
par(oma=c(1,2,1,0))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
hist(mean_BK_indCall,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,200000),xlim=c(-1.85,2))
mgp.axis(pos=-2,side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=axTicks(side=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"CNV call",line=-2,outer=T,cex=1.5)

hist(mean_BK_normcov,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,600000))
mgp.axis(pos=min(hist(mean_BK_normcov,breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=formatC(axTicks(side=1),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"Coverage",line=-2,outer=T,cex=1.5)

dev.off()


pdf("Hist_Cnmops_values_BukoKowa_paper_square_log.pdf",width=10,height=8,paper="special",pointsize=28)
par(oma=c(1,2,1,0))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
hist(mean_BK_indCall,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,200000),xlim=c(-1.85,2))
mgp.axis(pos=-2,side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=axTicks(side=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"CNV call",line=-2,outer=T,cex=1.5)

options(scipen=999)
hist(log10(mean_BK_normcov),breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,100000))
mgp.axis(pos=min(hist(log10(mean_BK_normcov),breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
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
HMA4_BK<-BK2[BK2$Scaffold==scaffold&BK2$Start>=(start-windowsize)&BK2$End<=(stop+windowsize),]
HMA4_KB<-KB2[KB2$Scaffold==scaffold&KB2$Start>=(start-windowsize)&KB2$End<=(stop+windowsize),]
HMA4_BK$Mid<-HMA4_BK$Start+((HMA4_BK$End-HMA4_BK$Start)/2)
HMA4_KB$Mid<-HMA4_KB$Start+((HMA4_BK$End-HMA4_BK$Start)/2)


require(Hmisc)
#HMA4
#scaffold_3	23475941 23483785	9845
BK_indCall2<-data.frame(Pos,BK_indCall)
KB_indCall2<-data.frame(Pos,KB_indCall)
BK_normcov2<-data.frame(Pos,BK_normcov)
KB_normcov2<-data.frame(Pos,KB_normcov)

HMA4_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
HMA4_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
HMA4_BK_indCall$Mid<-HMA4_BK_indCall$Start+((HMA4_BK_indCall$End-HMA4_BK_indCall$Start)/2)
HMA4_KB_indCall$Mid<-HMA4_KB_indCall$Start+((HMA4_BK_indCall$End-HMA4_BK_indCall$Start)/2)

HMA4_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
HMA4_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
HMA4_BK_normcov$Mid<-HMA4_BK_normcov$Start+((HMA4_BK_normcov$End-HMA4_BK_normcov$Start)/2)
HMA4_KB_normcov$Mid<-HMA4_KB_normcov$Start+((HMA4_BK_normcov$End-HMA4_BK_normcov$Start)/2)



require(Hmisc)
#HMA4
#scaffold_3	23475941 23483785	9845
pdf("Cnmops_Buko_HMA4_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(HMA4_BK_indCall[,4]~HMA4_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA4_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],HMA4_KB_indCall[,4:11]),max(HMA4_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],HMA4_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(HMA4_BK_indCall[,i]~HMA4_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA4_KB_indCall[,j]~HMA4_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA4_BK_normcov[,4]~HMA4_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],HMA4_KB_normcov[,4:11]),max(HMA4_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],HMA4_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(HMA4_BK_normcov[,i]~HMA4_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA4_KB_normcov[,j]~HMA4_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#MTP1
scaffold="scaffold_4"
start=22776960
stop=22778882
windowsize=50000

MTP1_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
MTP1_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
MTP1_BK_indCall$Mid<-MTP1_BK_indCall$Start+((MTP1_BK_indCall$End-MTP1_BK_indCall$Start)/2)
MTP1_KB_indCall$Mid<-MTP1_KB_indCall$Start+((MTP1_BK_indCall$End-MTP1_BK_indCall$Start)/2)

MTP1_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
MTP1_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
MTP1_BK_normcov$Mid<-MTP1_BK_normcov$Start+((MTP1_BK_normcov$End-MTP1_BK_normcov$Start)/2)
MTP1_KB_normcov$Mid<-MTP1_KB_normcov$Start+((MTP1_BK_normcov$End-MTP1_BK_normcov$Start)/2)

pdf("Cnmops_Buko_MTP1_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(MTP1_BK_indCall[,4]~MTP1_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(MTP1_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],MTP1_KB_indCall[,4:11]),max(MTP1_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],MTP1_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(MTP1_BK_indCall[,i]~MTP1_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(MTP1_KB_indCall[,j]~MTP1_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(MTP1_BK_normcov[,4]~MTP1_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(MTP1_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],MTP1_KB_normcov[,4:11]),max(MTP1_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],MTP1_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(MTP1_BK_normcov[,i]~MTP1_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(MTP1_KB_normcov[,j]~MTP1_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


#HMA3
scaffold="scaffold_7"
start=4960535
stop=4963880
windowsize=50000

HMA3_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
HMA3_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
HMA3_BK_indCall$Mid<-HMA3_BK_indCall$Start+((HMA3_BK_indCall$End-HMA3_BK_indCall$Start)/2)
HMA3_KB_indCall$Mid<-HMA3_KB_indCall$Start+((HMA3_BK_indCall$End-HMA3_BK_indCall$Start)/2)

HMA3_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
HMA3_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
HMA3_BK_normcov$Mid<-HMA3_BK_normcov$Start+((HMA3_BK_normcov$End-HMA3_BK_normcov$Start)/2)
HMA3_KB_normcov$Mid<-HMA3_KB_normcov$Start+((HMA3_BK_normcov$End-HMA3_BK_normcov$Start)/2)

pdf("Cnmops_Buko_HMA3_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(HMA3_BK_indCall[,4]~HMA3_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA3_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],HMA3_KB_indCall[,4:11]),max(HMA3_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],HMA3_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(HMA3_BK_indCall[,i]~HMA3_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA3_KB_indCall[,j]~HMA3_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA3_BK_normcov[,4]~HMA3_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA3_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],HMA3_KB_normcov[,4:11]),max(HMA3_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],HMA3_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(HMA3_BK_normcov[,i]~HMA3_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA3_KB_normcov[,j]~HMA3_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP6
scaffold="scaffold_4"
start=13949254
stop=13950804
windowsize=50000

ZIP6_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
ZIP6_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
ZIP6_BK_indCall$Mid<-ZIP6_BK_indCall$Start+((ZIP6_BK_indCall$End-ZIP6_BK_indCall$Start)/2)
ZIP6_KB_indCall$Mid<-ZIP6_KB_indCall$Start+((ZIP6_BK_indCall$End-ZIP6_BK_indCall$Start)/2)

ZIP6_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
ZIP6_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
ZIP6_BK_normcov$Mid<-ZIP6_BK_normcov$Start+((ZIP6_BK_normcov$End-ZIP6_BK_normcov$Start)/2)
ZIP6_KB_normcov$Mid<-ZIP6_KB_normcov$Start+((ZIP6_BK_normcov$End-ZIP6_BK_normcov$Start)/2)

pdf("Cnmops_Buko_ZIP6_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(ZIP6_BK_indCall[,4]~ZIP6_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP6_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],ZIP6_KB_indCall[,4:11]),max(ZIP6_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],ZIP6_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(ZIP6_BK_indCall[,i]~ZIP6_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP6_KB_indCall[,j]~ZIP6_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP6_BK_normcov[,4]~ZIP6_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP6_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],ZIP6_KB_normcov[,4:11]),max(ZIP6_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],ZIP6_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(ZIP6_BK_normcov[,i]~ZIP6_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP6_KB_normcov[,j]~ZIP6_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#NRAMP3
scaffold="scaffold_4"
start=1634727
stop=1637073
windowsize=50000

NRAMP3_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
NRAMP3_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
NRAMP3_BK_indCall$Mid<-NRAMP3_BK_indCall$Start+((NRAMP3_BK_indCall$End-NRAMP3_BK_indCall$Start)/2)
NRAMP3_KB_indCall$Mid<-NRAMP3_KB_indCall$Start+((NRAMP3_BK_indCall$End-NRAMP3_BK_indCall$Start)/2)

NRAMP3_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
NRAMP3_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
NRAMP3_BK_normcov$Mid<-NRAMP3_BK_normcov$Start+((NRAMP3_BK_normcov$End-NRAMP3_BK_normcov$Start)/2)
NRAMP3_KB_normcov$Mid<-NRAMP3_KB_normcov$Start+((NRAMP3_BK_normcov$End-NRAMP3_BK_normcov$Start)/2)

pdf("Cnmops_Buko_NRAMP3_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(NRAMP3_BK_indCall[,4]~NRAMP3_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(NRAMP3_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],NRAMP3_KB_indCall[,4:11]),max(NRAMP3_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],NRAMP3_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,2000,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(NRAMP3_BK_indCall[,i]~NRAMP3_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(NRAMP3_KB_indCall[,j]~NRAMP3_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(NRAMP3_BK_normcov[,4]~NRAMP3_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(NRAMP3_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],NRAMP3_KB_normcov[,4:11]),max(NRAMP3_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],NRAMP3_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,2000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(NRAMP3_BK_normcov[,i]~NRAMP3_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(NRAMP3_KB_normcov[,j]~NRAMP3_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP11
scaffold="scaffold_1"
start=29573151
stop=29574444
windowsize=50000

ZIP11_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
ZIP11_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
ZIP11_BK_indCall$Mid<-ZIP11_BK_indCall$Start+((ZIP11_BK_indCall$End-ZIP11_BK_indCall$Start)/2)
ZIP11_KB_indCall$Mid<-ZIP11_KB_indCall$Start+((ZIP11_BK_indCall$End-ZIP11_BK_indCall$Start)/2)

ZIP11_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
ZIP11_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
ZIP11_BK_normcov$Mid<-ZIP11_BK_normcov$Start+((ZIP11_BK_normcov$End-ZIP11_BK_normcov$Start)/2)
ZIP11_KB_normcov$Mid<-ZIP11_KB_normcov$Start+((ZIP11_BK_normcov$End-ZIP11_BK_normcov$Start)/2)

pdf("Cnmops_Buko_ZIP11_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(ZIP11_BK_indCall[,4]~ZIP11_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP11_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],ZIP11_KB_indCall[,4:11]),max(ZIP11_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],ZIP11_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(ZIP11_BK_indCall[,i]~ZIP11_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP11_KB_indCall[,j]~ZIP11_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP11_BK_normcov[,4]~ZIP11_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP11_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],ZIP11_KB_normcov[,4:11]),max(ZIP11_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],ZIP11_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(ZIP11_BK_normcov[,i]~ZIP11_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP11_KB_normcov[,j]~ZIP11_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#HMA2
scaffold="scaffold_7"
start=4971705
stop=4976723
windowsize=50000

HMA2_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
HMA2_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
HMA2_BK_indCall$Mid<-HMA2_BK_indCall$Start+((HMA2_BK_indCall$End-HMA2_BK_indCall$Start)/2)
HMA2_KB_indCall$Mid<-HMA2_KB_indCall$Start+((HMA2_BK_indCall$End-HMA2_BK_indCall$Start)/2)

HMA2_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
HMA2_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
HMA2_BK_normcov$Mid<-HMA2_BK_normcov$Start+((HMA2_BK_normcov$End-HMA2_BK_normcov$Start)/2)
HMA2_KB_normcov$Mid<-HMA2_KB_normcov$Start+((HMA2_BK_normcov$End-HMA2_BK_normcov$Start)/2)

pdf("Cnmops_Buko_HMA2_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(HMA2_BK_indCall[,4]~HMA2_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA2_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],HMA2_KB_indCall[,4:11]),max(HMA2_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],HMA2_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(HMA2_BK_indCall[,i]~HMA2_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA2_KB_indCall[,j]~HMA2_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA2_BK_normcov[,4]~HMA2_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA2_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],HMA2_KB_normcov[,4:11]),max(HMA2_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],HMA2_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(HMA2_BK_normcov[,i]~HMA2_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA2_KB_normcov[,j]~HMA2_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#IRT2
scaffold="scaffold_7"
start=10093506
stop=10094856
windowsize=50000

IRT2_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
IRT2_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
IRT2_BK_indCall$Mid<-IRT2_BK_indCall$Start+((IRT2_BK_indCall$End-IRT2_BK_indCall$Start)/2)
IRT2_KB_indCall$Mid<-IRT2_KB_indCall$Start+((IRT2_BK_indCall$End-IRT2_BK_indCall$Start)/2)

IRT2_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
IRT2_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
IRT2_BK_normcov$Mid<-IRT2_BK_normcov$Start+((IRT2_BK_normcov$End-IRT2_BK_normcov$Start)/2)
IRT2_KB_normcov$Mid<-IRT2_KB_normcov$Start+((IRT2_BK_normcov$End-IRT2_BK_normcov$Start)/2)

pdf("Cnmops_Buko_IRT2_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(IRT2_BK_indCall[,4]~IRT2_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(IRT2_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],IRT2_KB_indCall[,4:11]),max(IRT2_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],IRT2_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(IRT2_BK_indCall[,i]~IRT2_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(IRT2_KB_indCall[,j]~IRT2_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(IRT2_BK_normcov[,4]~IRT2_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(IRT2_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],IRT2_KB_normcov[,4:11]),max(IRT2_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],IRT2_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(IRT2_BK_normcov[,i]~IRT2_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(IRT2_KB_normcov[,j]~IRT2_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP8
scaffold="scaffold_8"
start=2517313
stop=2519529
windowsize=50000

ZIP8_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
ZIP8_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
ZIP8_BK_indCall$Mid<-ZIP8_BK_indCall$Start+((ZIP8_BK_indCall$End-ZIP8_BK_indCall$Start)/2)
ZIP8_KB_indCall$Mid<-ZIP8_KB_indCall$Start+((ZIP8_BK_indCall$End-ZIP8_BK_indCall$Start)/2)

ZIP8_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
ZIP8_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
ZIP8_BK_normcov$Mid<-ZIP8_BK_normcov$Start+((ZIP8_BK_normcov$End-ZIP8_BK_normcov$Start)/2)
ZIP8_KB_normcov$Mid<-ZIP8_KB_normcov$Start+((ZIP8_BK_normcov$End-ZIP8_BK_normcov$Start)/2)

pdf("Cnmops_Buko_ZIP8_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(ZIP8_BK_indCall[,4]~ZIP8_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP8_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],ZIP8_KB_indCall[,4:11]),max(ZIP8_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],ZIP8_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(ZIP8_BK_indCall[,i]~ZIP8_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP8_KB_indCall[,j]~ZIP8_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP8_BK_normcov[,4]~ZIP8_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP8_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],ZIP8_KB_normcov[,4:11]),max(ZIP8_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],ZIP8_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(ZIP8_BK_normcov[,i]~ZIP8_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP8_KB_normcov[,j]~ZIP8_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


#FRO5
scaffold="scaffold_6"
start=10296305
stop=10299581
windowsize=50000

FRO5_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
FRO5_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
FRO5_BK_indCall$Mid<-FRO5_BK_indCall$Start+((FRO5_BK_indCall$End-FRO5_BK_indCall$Start)/2)
FRO5_KB_indCall$Mid<-FRO5_KB_indCall$Start+((FRO5_BK_indCall$End-FRO5_BK_indCall$Start)/2)

FRO5_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
FRO5_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
FRO5_BK_normcov$Mid<-FRO5_BK_normcov$Start+((FRO5_BK_normcov$End-FRO5_BK_normcov$Start)/2)
FRO5_KB_normcov$Mid<-FRO5_KB_normcov$Start+((FRO5_BK_normcov$End-FRO5_BK_normcov$Start)/2)

pdf("Cnmops_Buko_FRO5_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(FRO5_BK_indCall[,4]~FRO5_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FRO5_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],FRO5_KB_indCall[,4:11]),max(FRO5_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],FRO5_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(FRO5_BK_indCall[,i]~FRO5_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRO5_KB_indCall[,j]~FRO5_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FRO5_BK_normcov[,4]~FRO5_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FRO5_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],FRO5_KB_normcov[,4:11]),max(FRO5_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],FRO5_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(FRO5_BK_normcov[,i]~FRO5_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRO5_KB_normcov[,j]~FRO5_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#FER1
scaffold="scaffold_6"
start=187631
stop=189383
windowsize=50000

FER1_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
FER1_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
FER1_BK_indCall$Mid<-FER1_BK_indCall$Start+((FER1_BK_indCall$End-FER1_BK_indCall$Start)/2)
FER1_KB_indCall$Mid<-FER1_KB_indCall$Start+((FER1_BK_indCall$End-FER1_BK_indCall$Start)/2)

FER1_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
FER1_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
FER1_BK_normcov$Mid<-FER1_BK_normcov$Start+((FER1_BK_normcov$End-FER1_BK_normcov$Start)/2)
FER1_KB_normcov$Mid<-FER1_KB_normcov$Start+((FER1_BK_normcov$End-FER1_BK_normcov$Start)/2)

pdf("Cnmops_Buko_FER1_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(FER1_BK_indCall[,4]~FER1_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FER1_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],FER1_KB_indCall[,4:11]),max(FER1_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],FER1_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(FER1_BK_indCall[,i]~FER1_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FER1_KB_indCall[,j]~FER1_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FER1_BK_normcov[,4]~FER1_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FER1_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],FER1_KB_normcov[,4:11]),max(FER1_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],FER1_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(FER1_BK_normcov[,i]~FER1_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FER1_KB_normcov[,j]~FER1_KB_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#FRD3
scaffold="scaffold_3"
start=3352004
stop=3357699
windowsize=50000

FRD3_BK_indCall<-BK_indCall2[BK_indCall2$Scaffold==scaffold&BK_indCall2$Start>=(start-windowsize)&BK_indCall2$End<=(stop+windowsize),]
FRD3_KB_indCall<-KB_indCall2[KB_indCall2$Scaffold==scaffold&KB_indCall2$Start>=(start-windowsize)&KB_indCall2$End<=(stop+windowsize),]
FRD3_BK_indCall$Mid<-FRD3_BK_indCall$Start+((FRD3_BK_indCall$End-FRD3_BK_indCall$Start)/2)
FRD3_KB_indCall$Mid<-FRD3_KB_indCall$Start+((FRD3_BK_indCall$End-FRD3_BK_indCall$Start)/2)

FRD3_BK_normcov<-BK_normcov2[BK_normcov2$Scaffold==scaffold&BK_normcov2$Start>=(start-windowsize)&BK_normcov2$End<=(stop+windowsize),]
FRD3_KB_normcov<-KB_normcov2[KB_normcov2$Scaffold==scaffold&KB_normcov2$Start>=(start-windowsize)&KB_normcov2$End<=(stop+windowsize),]
FRD3_BK_normcov$Mid<-FRD3_BK_normcov$Start+((FRD3_BK_normcov$End-FRD3_BK_normcov$Start)/2)
FRD3_KB_normcov$Mid<-FRD3_KB_normcov$Start+((FRD3_BK_normcov$End-FRD3_BK_normcov$Start)/2)

pdf("Cnmops_Buko_FRD3_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(FRD3_BK_indCall[,4]~FRD3_BK_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FRD3_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],FRD3_KB_indCall[,4:11]),max(FRD3_BK_indCall[,4:(ncol(HMA4_BK_indCall)-1)],FRD3_KB_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(FRD3_BK_indCall[,i]~FRD3_BK_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRD3_KB_indCall[,j]~FRD3_KB_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FRD3_BK_normcov[,4]~FRD3_BK_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FRD3_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],FRD3_KB_normcov[,4:11]),max(FRD3_BK_normcov[,4:(ncol(HMA4_BK_indCall)-1)],FRD3_KB_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_BK_indCall)-1))
	{lines(FRD3_BK_normcov[,i]~FRD3_BK_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRD3_KB_normcov[,j]~FRD3_KB_normcov$Mid,col="black",lwd=2.5)}
box()s
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Buko"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


