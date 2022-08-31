PK_BL_topo3<-read.delim("PiekhalleriKowalyr_w100_branchlengths.topo3.txt",sep=" ",header=F)
PK_BL_topo3<-PK_BL_topo3[!PK_BL_topo3$V1=="0",]
PK_BL_topo2<-read.delim("PiekhalleriKowalyr_w100_branchlengths.topo2.txt",sep=" ",header=F)
PK_BL_topo2<-PK_BL_topo2[!PK_BL_topo2$V1=="0",]
PK_BL_topo1<-read.delim("PiekhalleriKowalyr_w100_branchlengths.topo1.txt",sep=" ",header=F)
PK_BL_topo1<-PK_BL_topo1[!PK_BL_topo1$V1=="0",]

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
Modes(PK_BL_topo1$V5)#halleri from Kowaar
#1.01
Modes(PK_BL_topo2$V5)#Piekar from Kowaar
#0.88
Modes(PK_BL_topo3$V5)#Piekar from halleri
#0.74 0.81

mean(PK_BL_topo1$V5)
#1.076379
mean(PK_BL_topo2$V5)
#0.8819987
mean(PK_BL_topo3$V5)
#0.8255016

median(PK_BL_topo1$V5)
#1.02
median(PK_BL_topo2$V5)
#0.87
median(PK_BL_topo3$V5)
#0.8

sd(PK_BL_topo1$V5)
#0.3882837
sd(PK_BL_topo2$V5)
#0.2957407
sd(PK_BL_topo3$V5)
#0.361842




MK_BL_topo3<-read.delim("MiashalleriKowaarenosalyr_branchlengths_topo3_w100.txt",sep=" ",header=F)
MK_BL_topo3<-MK_BL_topo3[!MK_BL_topo3$V1=="0",]
MK_BL_topo2<-read.delim("MiashalleriKowaarenosalyr_branchlengths_topo2_w100.txt",sep=" ",header=F)
MK_BL_topo2<-MK_BL_topo2[!MK_BL_topo2$V1=="0",]
MK_BL_topo1<-read.delim("MiashalleriKowaarenosalyr_branchlengths_topo1_w100.txt",sep=" ",header=F)
MK_BL_topo1<-MK_BL_topo1[!MK_BL_topo1$V1=="0",]

Modes(MK_BL_topo1$V5)#halleri from Kowaar
#1.08
Modes(MK_BL_topo2$V5)#Piekar from Kowaar
#0.8
Modes(MK_BL_topo3$V5)#Piekar from halleri
#0.83

mean(MK_BL_topo1$V5)
#1.077771
mean(MK_BL_topo2$V5)
#0.8761776
mean(MK_BL_topo3$V5)
#0.8606867

median(MK_BL_topo1$V5)
#1.03
median(MK_BL_topo2$V5)
#0.86
median(MK_BL_topo3$V5)
#0.83

sd(MK_BL_topo1$V5)
#0.3856965
sd(MK_BL_topo2$V5)
#0.2922455
sd(MK_BL_topo3$V5)
#0.3542425

BK_BL_topo3<-read.delim("BukohalleriKowalyr_w100_branchlengths.topo3.txt",sep=" ",header=F)
BK_BL_topo3<-BK_BL_topo3[!BK_BL_topo3$V1=="0",]
BK_BL_topo2<-read.delim("BukohalleriKowalyr_w100_branchlengths.topo2.txt",sep=" ",header=F)
BK_BL_topo2<-BK_BL_topo2[!BK_BL_topo2$V1=="0",]
BK_BL_topo1<-read.delim("BukohalleriKowalyr_w100_branchlengths.topo1.txt",sep=" ",header=F)
BK_BL_topo1<-BK_BL_topo1[!BK_BL_topo1$V1=="0",]

Modes(BK_BL_topo1$V5)#halleri from Kowaar
#0.9
Modes(BK_BL_topo2$V5)#Bukoar from Kowaar
#0.88
Modes(BK_BL_topo3$V5)#Bukoar from halleri
#0.77

mean(BK_BL_topo1$V5)
#1.08572
mean(BK_BL_topo2$V5)
#0.9111552
mean(BK_BL_topo3$V5)
#0.8919999

median(BK_BL_topo1$V5)
#1.04
median(BK_BL_topo2$V5)
#0.89
median(BK_BL_topo3$V5)
#0.86

sd(BK_BL_topo1$V5)
#0.393644
sd(BK_BL_topo2$V5)
#0.3050194
sd(BK_BL_topo3$V5)
#0.3741623


KK_BL_topo3<-read.delim("KatohalleriKowalyr_w100_branchlengths.topo3.txt",sep=" ",header=F)
KK_BL_topo3<-KK_BL_topo3[!KK_BL_topo3$V1=="0",]
KK_BL_topo2<-read.delim("KatohalleriKowalyr_w100_branchlengths.topo2.txt",sep=" ",header=F)
KK_BL_topo2<-KK_BL_topo2[!KK_BL_topo2$V1=="0",]
KK_BL_topo1<-read.delim("KatohalleriKowalyr_w100_branchlengths.topo1.txt",sep=" ",header=F)
KK_BL_topo1<-KK_BL_topo1[!KK_BL_topo1$V1=="0",]

Modes(KK_BL_topo1$V5)#halleri from Kowaar
#0.93
Modes(KK_BL_topo2$V5)#Katoar from Kowaar
#0.85
Modes(KK_BL_topo3$V5)#Katoar from halleri
#0.82

mean(KK_BL_topo1$V5)
#1.080627
mean(KK_BL_topo2$V5)
#0.8946249
mean(KK_BL_topo3$V5)
#0.8651893

median(KK_BL_topo1$V5)
#1.03
median(KK_BL_topo2$V5)
#0.88
median(KK_BL_topo3$V5)
#0.84

sd(KK_BL_topo1$V5)
#0.3893761
sd(KK_BL_topo2$V5)
#0.301388
sd(KK_BL_topo3$V5)
#0.3632154

pdf("Twisst_branch lengths_w100_paper.pdf",width=15,height=8,paper="special",pointsize=10)
par(mfrow=c(2,2))
par(mgp=c(2.5,0.75,0))
par(oma=c(2,2,2,2))
par(lwd=2.5)
par(mar=c(3,5,1,1))
boxplot(MK_BL_topo1$V5,MK_BL_topo2$V5,MK_BL_topo3$V5,xlim=c(0.6,3.4),xlab="",ylab="Distance",col="white",border=c("blue","black","red"),names=c("Topology 1","Topology 2","Topology 3"),cex.axis=1.3,cex.lab=1.8,las=1,lwd=2.5)
boxplot(PK_BL_topo1$V5,PK_BL_topo2$V5,PK_BL_topo3$V5,xlim=c(0.6,3.4),xlab="",ylab="",col="white",border=c("blue","black","red"),names=c("Topology 1","Topology 2","Topology 3"),cex.axis=1.3,cex.lab=1.8,las=1,lwd=2.5)
boxplot(BK_BL_topo1$V5,BK_BL_topo2$V5,BK_BL_topo3$V5,xlim=c(0.6,3.4),xlab="",ylab="Distance",col="white",border=c("blue","black","red"),names=c("Topology 1","Topology 2","Topology 3"),cex.axis=1.3,cex.lab=1.8,las=1,lwd=2.5)
boxplot(KK_BL_topo1$V5,KK_BL_topo2$V5,KK_BL_topo3$V5,xlim=c(0.6,3.4),xlab="",ylab="",col="white",border=c("blue","black","red"),names=c("Topology 1","Topology 2","Topology 3"),cex.axis=1.3,cex.lab=1.8,las=1,lwd=2.5)

dev.off()




