require(openxlsx)
KaKo_indCall<-read.xlsx("IndividualCall_data_cnmops_KaKoar.xlsx",1)
KaKo_normcov<-read.xlsx("Normalized_data_cnmops_KaKoar.xlsx",1)

KoKa_indCall<-read.xlsx("IndividualCall_data_cnmops_KKar.xlsx",1)
KoKa_normcov<-read.xlsx("Normalized_data_cnmops_KKar.xlsx",1)

mean_KoKa_indCall<-apply(KoKa_indCall,1,mean)
mean_KoKa_normcov<-apply(KoKa_normcov,1,mean)
mean_KaKo_indCall<-apply(KaKo_indCall,1,mean)
mean_KaKo_normcov<-apply(KaKo_normcov,1,mean)

KaKo<-data.frame(mean_KaKo_indCall,mean_KaKo_normcov)
KoKa<-data.frame(mean_KoKa_indCall,mean_KoKa_normcov)

Pos<-read.xlsx("Pos__data_cnmops_MKar.xlsx",1)
require(tidyr)
Pos<-separate(Pos,1,sep="_",into=c("NA","Scaff","Start","End"),remove=T)
Pos$Scaffold<-paste(Pos[,1],Pos[,2],sep="_")
Pos<-Pos[,c(5,3,4)]
Pos[,2]<-as.numeric(Pos[,2])
Pos[,3]<-as.numeric(Pos[,3])

KaKo2<-data.frame(Pos,KaKo)
KoKa2<-data.frame(Pos,KoKa)

require(dichromat)
colfunc <- colorRampPalette(c("black","red"))
colours_10<-colfunc(10)


pdf("Cnmops_values.pdf",width=8, height=12,paper="special")
par(mar=c(6,6,0,1)+0.1)
par(mfrow=c(3,1))
plot(mean_KoKa_locass,mean_KoKa_indCall)
plot(mean_KoKa_normcov,mean_KoKa_indCall)
plot(mean_KoKa_posteriorprobs_CN1,mean_KoKa_indCall,col=colours_10[1])
for (i in 2:8)
	{points(get(paste0("mean_KoKa_posteriorprobs_CN",i)),mean_KoKa_indCall,col=colours_10[i])}
dev.off()

pdf("Cnmops_values_hist.pdf",width=12, height=12,paper="special")
par(mar=c(6,6,2,1)+0.1)
par(mfrow=c(2,2))
hist(mean_KoKa_indCall)
hist(mean_KoKa_locass)
hist(mean_KoKa_normcov)
hist(mean_KoKa_posteriorprobs_CN1,col=colours_10[1])
for (i in 2:8)
	{hist(get(paste0("mean_KoKa_posteriorprobs_CN",i)),col=colours_10[i],add=T)}
dev.off()
options(scipen=999)
require(Hmisc)
pdf("Hist_Cnmops_values_KatoKowa_paper_square.pdf",width=10,height=8,paper="special",pointsize=25)
par(oma=c(1,2,1,0))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
hist(mean_KaKo_indCall,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,200000),xlim=c(-1.85,2))
mgp.axis(pos=-2,side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=axTicks(side=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"CNV call",line=-2,outer=T,cex=1.5)

hist(mean_KaKo_normcov,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,600000))
mgp.axis(pos=min(hist(mean_KaKo_normcov,breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=formatC(axTicks(side=1),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"Coverage",line=-2,outer=T,cex=1.5)

dev.off()


pdf("Hist_Cnmops_values_KatoKowa_paper_square_log.pdf",width=10,height=8,paper="special",pointsize=28)
par(oma=c(1,2,1,0))
par(mgp=c(2.5,0.75,0))
par(mar=c(4,4,1,1))
hist(mean_KaKo_indCall,breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,200000),xlim=c(-1.85,2))
mgp.axis(pos=-2,side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mgp.axis(pos=0,side=1,at=axTicks(side=1),labels=axTicks(side=1),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
mtext(side=2,"Frequency",line=1,outer=T,cex=1.5)
mtext(side=1,"CNV call",line=-2,outer=T,cex=1.5)

options(scipen=999)
hist(log10(mean_KaKo_normcov),breaks=100,main="",cex.lab=1.5,cex.main=1.5,font.main=1.5,xlab="",las=1,yaxt="n",ylab="",col="black",xaxt="n",ylim=c(0,100000))
mgp.axis(pos=min(hist(log10(mean_KaKo_normcov),breaks=100,plot=F)$breaks),side=2,at=axTicks(side=2),labels=formatC(axTicks(side=2),big.mark=",",format="g",decimal.mark = ".",digits=6),cex=1,cex.lab=1.2,cex.axis=1.3,las=1)
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
HMA4_KaKo<-KaKo2[KaKo2$Scaffold==scaffold&KaKo2$Start>=(start-windowsize)&KaKo2$End<=(stop+windowsize),]
HMA4_KoKa<-KoKa2[KoKa2$Scaffold==scaffold&KoKa2$Start>=(start-windowsize)&KoKa2$End<=(stop+windowsize),]
HMA4_KaKo$Mid<-HMA4_KaKo$Start+((HMA4_KaKo$End-HMA4_KaKo$Start)/2)
HMA4_KoKa$Mid<-HMA4_KoKa$Start+((HMA4_KaKo$End-HMA4_KaKo$Start)/2)


require(Hmisc)
#HMA4
#scaffold_3	23475941 23483785	9845
KaKo_indCall2<-data.frame(Pos,KaKo_indCall)
KoKa_indCall2<-data.frame(Pos,KoKa_indCall)
KaKo_normcov2<-data.frame(Pos,KaKo_normcov)
KoKa_normcov2<-data.frame(Pos,KoKa_normcov)

HMA4_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
HMA4_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
HMA4_KaKo_indCall$Mid<-HMA4_KaKo_indCall$Start+((HMA4_KaKo_indCall$End-HMA4_KaKo_indCall$Start)/2)
HMA4_KoKa_indCall$Mid<-HMA4_KoKa_indCall$Start+((HMA4_KaKo_indCall$End-HMA4_KaKo_indCall$Start)/2)

HMA4_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
HMA4_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
HMA4_KaKo_normcov$Mid<-HMA4_KaKo_normcov$Start+((HMA4_KaKo_normcov$End-HMA4_KaKo_normcov$Start)/2)
HMA4_KoKa_normcov$Mid<-HMA4_KoKa_normcov$Start+((HMA4_KaKo_normcov$End-HMA4_KaKo_normcov$Start)/2)



require(Hmisc)
#HMA4
#scaffold_3	23475941 23483785	9845
pdf("Cnmops_Kato_HMA4_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(HMA4_KaKo_indCall[,4]~HMA4_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA4_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA4_KoKa_indCall[,4:11]),max(HMA4_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA4_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(HMA4_KaKo_indCall[,i]~HMA4_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA4_KoKa_indCall[,j]~HMA4_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA4_KaKo_normcov[,4]~HMA4_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA4_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA4_KoKa_normcov[,4:11]),max(HMA4_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA4_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(HMA4_KaKo_normcov[,i]~HMA4_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA4_KoKa_normcov[,j]~HMA4_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#MTP1
scaffold="scaffold_4"
start=22776960
stop=22778882
windowsize=50000

MTP1_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
MTP1_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
MTP1_KaKo_indCall$Mid<-MTP1_KaKo_indCall$Start+((MTP1_KaKo_indCall$End-MTP1_KaKo_indCall$Start)/2)
MTP1_KoKa_indCall$Mid<-MTP1_KoKa_indCall$Start+((MTP1_KaKo_indCall$End-MTP1_KaKo_indCall$Start)/2)

MTP1_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
MTP1_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
MTP1_KaKo_normcov$Mid<-MTP1_KaKo_normcov$Start+((MTP1_KaKo_normcov$End-MTP1_KaKo_normcov$Start)/2)
MTP1_KoKa_normcov$Mid<-MTP1_KoKa_normcov$Start+((MTP1_KaKo_normcov$End-MTP1_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_MTP1_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(MTP1_KaKo_indCall[,4]~MTP1_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(MTP1_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],MTP1_KoKa_indCall[,4:11]),max(MTP1_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],MTP1_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(MTP1_KaKo_indCall[,i]~MTP1_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(MTP1_KoKa_indCall[,j]~MTP1_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(MTP1_KaKo_normcov[,4]~MTP1_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(MTP1_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],MTP1_KoKa_normcov[,4:11]),max(MTP1_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],MTP1_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(MTP1_KaKo_normcov[,i]~MTP1_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(MTP1_KoKa_normcov[,j]~MTP1_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


#HMA3
scaffold="scaffold_7"
start=4960535
stop=4963880
windowsize=50000

HMA3_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
HMA3_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
HMA3_KaKo_indCall$Mid<-HMA3_KaKo_indCall$Start+((HMA3_KaKo_indCall$End-HMA3_KaKo_indCall$Start)/2)
HMA3_KoKa_indCall$Mid<-HMA3_KoKa_indCall$Start+((HMA3_KaKo_indCall$End-HMA3_KaKo_indCall$Start)/2)

HMA3_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
HMA3_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
HMA3_KaKo_normcov$Mid<-HMA3_KaKo_normcov$Start+((HMA3_KaKo_normcov$End-HMA3_KaKo_normcov$Start)/2)
HMA3_KoKa_normcov$Mid<-HMA3_KoKa_normcov$Start+((HMA3_KaKo_normcov$End-HMA3_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_HMA3_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(HMA3_KaKo_indCall[,4]~HMA3_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA3_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA3_KoKa_indCall[,4:11]),max(HMA3_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA3_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(HMA3_KaKo_indCall[,i]~HMA3_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA3_KoKa_indCall[,j]~HMA3_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA3_KaKo_normcov[,4]~HMA3_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA3_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA3_KoKa_normcov[,4:11]),max(HMA3_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA3_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(HMA3_KaKo_normcov[,i]~HMA3_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA3_KoKa_normcov[,j]~HMA3_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP6
scaffold="scaffold_4"
start=13949254
stop=13950804
windowsize=50000

ZIP6_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
ZIP6_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
ZIP6_KaKo_indCall$Mid<-ZIP6_KaKo_indCall$Start+((ZIP6_KaKo_indCall$End-ZIP6_KaKo_indCall$Start)/2)
ZIP6_KoKa_indCall$Mid<-ZIP6_KoKa_indCall$Start+((ZIP6_KaKo_indCall$End-ZIP6_KaKo_indCall$Start)/2)

ZIP6_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
ZIP6_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
ZIP6_KaKo_normcov$Mid<-ZIP6_KaKo_normcov$Start+((ZIP6_KaKo_normcov$End-ZIP6_KaKo_normcov$Start)/2)
ZIP6_KoKa_normcov$Mid<-ZIP6_KoKa_normcov$Start+((ZIP6_KaKo_normcov$End-ZIP6_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_ZIP6_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(ZIP6_KaKo_indCall[,4]~ZIP6_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP6_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP6_KoKa_indCall[,4:11]),max(ZIP6_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP6_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(ZIP6_KaKo_indCall[,i]~ZIP6_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP6_KoKa_indCall[,j]~ZIP6_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP6_KaKo_normcov[,4]~ZIP6_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP6_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP6_KoKa_normcov[,4:11]),max(ZIP6_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP6_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(ZIP6_KaKo_normcov[,i]~ZIP6_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP6_KoKa_normcov[,j]~ZIP6_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#NRAMP3
scaffold="scaffold_4"
start=1634727
stop=1637073
windowsize=50000

NRAMP3_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
NRAMP3_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
NRAMP3_KaKo_indCall$Mid<-NRAMP3_KaKo_indCall$Start+((NRAMP3_KaKo_indCall$End-NRAMP3_KaKo_indCall$Start)/2)
NRAMP3_KoKa_indCall$Mid<-NRAMP3_KoKa_indCall$Start+((NRAMP3_KaKo_indCall$End-NRAMP3_KaKo_indCall$Start)/2)

NRAMP3_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
NRAMP3_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
NRAMP3_KaKo_normcov$Mid<-NRAMP3_KaKo_normcov$Start+((NRAMP3_KaKo_normcov$End-NRAMP3_KaKo_normcov$Start)/2)
NRAMP3_KoKa_normcov$Mid<-NRAMP3_KoKa_normcov$Start+((NRAMP3_KaKo_normcov$End-NRAMP3_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_NRAMP3_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(NRAMP3_KaKo_indCall[,4]~NRAMP3_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(NRAMP3_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],NRAMP3_KoKa_indCall[,4:11]),max(NRAMP3_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],NRAMP3_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(NRAMP3_KaKo_indCall[,i]~NRAMP3_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(NRAMP3_KoKa_indCall[,j]~NRAMP3_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(NRAMP3_KaKo_normcov[,4]~NRAMP3_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(NRAMP3_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],NRAMP3_KoKa_normcov[,4:11]),max(NRAMP3_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],NRAMP3_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,2000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(NRAMP3_KaKo_normcov[,i]~NRAMP3_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(NRAMP3_KoKa_normcov[,j]~NRAMP3_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP11
scaffold="scaffold_1"
start=29573151
stop=29574444
windowsize=50000

ZIP11_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
ZIP11_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
ZIP11_KaKo_indCall$Mid<-ZIP11_KaKo_indCall$Start+((ZIP11_KaKo_indCall$End-ZIP11_KaKo_indCall$Start)/2)
ZIP11_KoKa_indCall$Mid<-ZIP11_KoKa_indCall$Start+((ZIP11_KaKo_indCall$End-ZIP11_KaKo_indCall$Start)/2)

ZIP11_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
ZIP11_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
ZIP11_KaKo_normcov$Mid<-ZIP11_KaKo_normcov$Start+((ZIP11_KaKo_normcov$End-ZIP11_KaKo_normcov$Start)/2)
ZIP11_KoKa_normcov$Mid<-ZIP11_KoKa_normcov$Start+((ZIP11_KaKo_normcov$End-ZIP11_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_ZIP11_paper_inds2.pdf",width=5,height=3,paper="special",pointsize=10)
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

plot(ZIP11_KaKo_indCall[,4]~ZIP11_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP11_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP11_KoKa_indCall[,4:11]),max(ZIP11_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP11_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(ZIP11_KaKo_indCall[,i]~ZIP11_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP11_KoKa_indCall[,j]~ZIP11_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP11_KaKo_normcov[,4]~ZIP11_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP11_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP11_KoKa_normcov[,4:11]),max(ZIP11_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP11_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(ZIP11_KaKo_normcov[,i]~ZIP11_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP11_KoKa_normcov[,j]~ZIP11_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#HMA2
scaffold="scaffold_7"
start=4971705
stop=4976723
windowsize=50000

HMA2_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
HMA2_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
HMA2_KaKo_indCall$Mid<-HMA2_KaKo_indCall$Start+((HMA2_KaKo_indCall$End-HMA2_KaKo_indCall$Start)/2)
HMA2_KoKa_indCall$Mid<-HMA2_KoKa_indCall$Start+((HMA2_KaKo_indCall$End-HMA2_KaKo_indCall$Start)/2)

HMA2_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
HMA2_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
HMA2_KaKo_normcov$Mid<-HMA2_KaKo_normcov$Start+((HMA2_KaKo_normcov$End-HMA2_KaKo_normcov$Start)/2)
HMA2_KoKa_normcov$Mid<-HMA2_KoKa_normcov$Start+((HMA2_KaKo_normcov$End-HMA2_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_HMA2_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(HMA2_KaKo_indCall[,4]~HMA2_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(HMA2_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA2_KoKa_indCall[,4:11]),max(HMA2_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA2_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(HMA2_KaKo_indCall[,i]~HMA2_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA2_KoKa_indCall[,j]~HMA2_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(HMA2_KaKo_normcov[,4]~HMA2_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(HMA2_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA2_KoKa_normcov[,4:11]),max(HMA2_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],HMA2_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(HMA2_KaKo_normcov[,i]~HMA2_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(HMA2_KoKa_normcov[,j]~HMA2_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#IRT2
scaffold="scaffold_7"
start=10093506
stop=10094856
windowsize=50000

IRT2_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
IRT2_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
IRT2_KaKo_indCall$Mid<-IRT2_KaKo_indCall$Start+((IRT2_KaKo_indCall$End-IRT2_KaKo_indCall$Start)/2)
IRT2_KoKa_indCall$Mid<-IRT2_KoKa_indCall$Start+((IRT2_KaKo_indCall$End-IRT2_KaKo_indCall$Start)/2)

IRT2_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
IRT2_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
IRT2_KaKo_normcov$Mid<-IRT2_KaKo_normcov$Start+((IRT2_KaKo_normcov$End-IRT2_KaKo_normcov$Start)/2)
IRT2_KoKa_normcov$Mid<-IRT2_KoKa_normcov$Start+((IRT2_KaKo_normcov$End-IRT2_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_IRT2_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(IRT2_KaKo_indCall[,4]~IRT2_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(IRT2_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],IRT2_KoKa_indCall[,4:11]),max(IRT2_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],IRT2_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(IRT2_KaKo_indCall[,i]~IRT2_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(IRT2_KoKa_indCall[,j]~IRT2_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(IRT2_KaKo_normcov[,4]~IRT2_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(IRT2_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],IRT2_KoKa_normcov[,4:11]),max(IRT2_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],IRT2_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(IRT2_KaKo_normcov[,i]~IRT2_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(IRT2_KoKa_normcov[,j]~IRT2_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#ZIP8
scaffold="scaffold_8"
start=2517313
stop=2519529
windowsize=50000

ZIP8_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
ZIP8_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
ZIP8_KaKo_indCall$Mid<-ZIP8_KaKo_indCall$Start+((ZIP8_KaKo_indCall$End-ZIP8_KaKo_indCall$Start)/2)
ZIP8_KoKa_indCall$Mid<-ZIP8_KoKa_indCall$Start+((ZIP8_KaKo_indCall$End-ZIP8_KaKo_indCall$Start)/2)

ZIP8_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
ZIP8_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
ZIP8_KaKo_normcov$Mid<-ZIP8_KaKo_normcov$Start+((ZIP8_KaKo_normcov$End-ZIP8_KaKo_normcov$Start)/2)
ZIP8_KoKa_normcov$Mid<-ZIP8_KoKa_normcov$Start+((ZIP8_KaKo_normcov$End-ZIP8_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_ZIP8_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(ZIP8_KaKo_indCall[,4]~ZIP8_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(ZIP8_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP8_KoKa_indCall[,4:11]),max(ZIP8_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP8_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(ZIP8_KaKo_indCall[,i]~ZIP8_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP8_KoKa_indCall[,j]~ZIP8_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(ZIP8_KaKo_normcov[,4]~ZIP8_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(ZIP8_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP8_KoKa_normcov[,4:11]),max(ZIP8_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],ZIP8_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(ZIP8_KaKo_normcov[,i]~ZIP8_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(ZIP8_KoKa_normcov[,j]~ZIP8_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


#FRO5
scaffold="scaffold_6"
start=10296305
stop=10299581
windowsize=50000

FRO5_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
FRO5_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
FRO5_KaKo_indCall$Mid<-FRO5_KaKo_indCall$Start+((FRO5_KaKo_indCall$End-FRO5_KaKo_indCall$Start)/2)
FRO5_KoKa_indCall$Mid<-FRO5_KoKa_indCall$Start+((FRO5_KaKo_indCall$End-FRO5_KaKo_indCall$Start)/2)

FRO5_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
FRO5_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
FRO5_KaKo_normcov$Mid<-FRO5_KaKo_normcov$Start+((FRO5_KaKo_normcov$End-FRO5_KaKo_normcov$Start)/2)
FRO5_KoKa_normcov$Mid<-FRO5_KoKa_normcov$Start+((FRO5_KaKo_normcov$End-FRO5_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_FRO5_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(FRO5_KaKo_indCall[,4]~FRO5_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FRO5_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],FRO5_KoKa_indCall[,4:11]),max(FRO5_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],FRO5_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(FRO5_KaKo_indCall[,i]~FRO5_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRO5_KoKa_indCall[,j]~FRO5_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FRO5_KaKo_normcov[,4]~FRO5_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FRO5_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],FRO5_KoKa_normcov[,4:11]),max(FRO5_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],FRO5_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(FRO5_KaKo_normcov[,i]~FRO5_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRO5_KoKa_normcov[,j]~FRO5_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#FER1
scaffold="scaffold_6"
start=187631
stop=189383
windowsize=50000

FER1_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
FER1_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
FER1_KaKo_indCall$Mid<-FER1_KaKo_indCall$Start+((FER1_KaKo_indCall$End-FER1_KaKo_indCall$Start)/2)
FER1_KoKa_indCall$Mid<-FER1_KoKa_indCall$Start+((FER1_KaKo_indCall$End-FER1_KaKo_indCall$Start)/2)

FER1_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
FER1_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
FER1_KaKo_normcov$Mid<-FER1_KaKo_normcov$Start+((FER1_KaKo_normcov$End-FER1_KaKo_normcov$Start)/2)
FER1_KoKa_normcov$Mid<-FER1_KoKa_normcov$Start+((FER1_KaKo_normcov$End-FER1_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_FER1_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(FER1_KaKo_indCall[,4]~FER1_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FER1_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],FER1_KoKa_indCall[,4:11]),max(FER1_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],FER1_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(FER1_KaKo_indCall[,i]~FER1_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FER1_KoKa_indCall[,j]~FER1_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FER1_KaKo_normcov[,4]~FER1_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FER1_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],FER1_KoKa_normcov[,4:11]),max(FER1_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],FER1_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(FER1_KaKo_normcov[,i]~FER1_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FER1_KoKa_normcov[,j]~FER1_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()

#FRD3
scaffold="scaffold_3"
start=3352004
stop=3357699
windowsize=50000

FRD3_KaKo_indCall<-KaKo_indCall2[KaKo_indCall2$Scaffold==scaffold&KaKo_indCall2$Start>=(start-windowsize)&KaKo_indCall2$End<=(stop+windowsize),]
FRD3_KoKa_indCall<-KoKa_indCall2[KoKa_indCall2$Scaffold==scaffold&KoKa_indCall2$Start>=(start-windowsize)&KoKa_indCall2$End<=(stop+windowsize),]
FRD3_KaKo_indCall$Mid<-FRD3_KaKo_indCall$Start+((FRD3_KaKo_indCall$End-FRD3_KaKo_indCall$Start)/2)
FRD3_KoKa_indCall$Mid<-FRD3_KoKa_indCall$Start+((FRD3_KaKo_indCall$End-FRD3_KaKo_indCall$Start)/2)

FRD3_KaKo_normcov<-KaKo_normcov2[KaKo_normcov2$Scaffold==scaffold&KaKo_normcov2$Start>=(start-windowsize)&KaKo_normcov2$End<=(stop+windowsize),]
FRD3_KoKa_normcov<-KoKa_normcov2[KoKa_normcov2$Scaffold==scaffold&KoKa_normcov2$Start>=(start-windowsize)&KoKa_normcov2$End<=(stop+windowsize),]
FRD3_KaKo_normcov$Mid<-FRD3_KaKo_normcov$Start+((FRD3_KaKo_normcov$End-FRD3_KaKo_normcov$Start)/2)
FRD3_KoKa_normcov$Mid<-FRD3_KoKa_normcov$Start+((FRD3_KaKo_normcov$End-FRD3_KaKo_normcov$Start)/2)

pdf("Cnmops_Kato_FRD3_paper_inds2.pdf",width=13,height=7,paper="special",pointsize=20)
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

plot(FRD3_KaKo_indCall[,4]~FRD3_KaKo_indCall$Mid,col="red",type="l",xlab="",ylab="",lwd=2.5,xaxt="n",cex.axis=1.3,las=1,cex.lab=1.5,ylim=c(min(FRD3_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],FRD3_KoKa_indCall[,4:11]),max(FRD3_KaKo_indCall[,4:(ncol(HMA4_KaKo_indCall)-1)],FRD3_KoKa_indCall[,4:11])))
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
rect(start,-100,stop,100,col="lightgrey",border = NA)	
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(FRD3_KaKo_indCall[,i]~FRD3_KaKo_indCall$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRD3_KoKa_indCall[,j]~FRD3_KoKa_indCall$Mid,col="black",lwd=2.5)}
mtext(side=2,"CNV call",line=3,xpd=T,cex=1.5,adj=0.3)
box()

#par(mar=c(2,6,1,1)+0.1)
plot(FRD3_KaKo_normcov[,4]~FRD3_KaKo_normcov$Mid,col="red",type="l",xlab="S",ylab="",lwd=2.5,xaxt="n",cex.axis=1.2,las=1,cex.lab=1.5,ylim=c(min(FRD3_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],FRD3_KoKa_normcov[,4:11]),max(FRD3_KaKo_normcov[,4:(ncol(HMA4_KaKo_indCall)-1)],FRD3_KoKa_normcov[,4:11])),xpd=T)
rect(start,-100,stop,1000,col="lightgrey",border = NA)	
mgp.axis(side=1,at=axTicks(side=1),labels=NA,cex=1,cex.lab=1.2,cex.axis=1.3,xpd=T)
for (i in 4:(ncol(HMA4_KaKo_indCall)-1))
	{lines(FRD3_KaKo_normcov[,i]~FRD3_KaKo_normcov$Mid,col="red",lwd=2.5)}
for (j in 4:11)
	{lines(FRD3_KoKa_normcov[,j]~FRD3_KoKa_normcov$Mid,col="black",lwd=2.5)}
box()
mtext(side=2,"Coverage",line=3,xpd=T,cex=1.5,adj=0.3)
#mtext(side=1,"Scaffold position (Mb)",line=2,outer=T,xpd=T,cex=1.5)
#legend("bottomleft",lwd=3,col=c("black","red"),legend=c("Kowa","Kato"),bty="n",cex=1.3,horiz=T,inset=c(-0.1,-0.4),xpd=TRUE)

dev.off()


