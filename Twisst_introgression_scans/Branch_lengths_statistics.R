MK_BL_topo3<-read.delim("MiashalleriKowaarenosalyr_branchlengths_topo3_w100.txt",sep=" ",header=F)
MK_BL_topo3<-MK_BL_topo3[!MK_BL_topo3$V1=="0",]
MK_BL_topo2<-read.delim("MiashalleriKowaarenosalyr_branchlengths_topo2_w100.txt",sep=" ",header=F)
MK_BL_topo2<-MK_BL_topo2[!MK_BL_topo2$V1=="0",]
MK_BL_topo1<-read.delim("MiashalleriKowaarenosalyr_branchlengths_topo1_w100.txt",sep=" ",header=F)
MK_BL_topo1<-MK_BL_topo1[!MK_BL_topo1$V1=="0",]

PK_BL_topo3<-read.delim("PiekhalleriKowalyr_w100_branchlengths.topo3.txt",sep=" ",header=F)
PK_BL_topo3<-PK_BL_topo3[!PK_BL_topo3$V1=="0",]
PK_BL_topo2<-read.delim("PiekhalleriKowalyr_w100_branchlengths.topo2.txt",sep=" ",header=F)
PK_BL_topo2<-PK_BL_topo2[!PK_BL_topo2$V1=="0",]
PK_BL_topo1<-read.delim("PiekhalleriKowalyr_w100_branchlengths.topo1.txt",sep=" ",header=F)
PK_BL_topo1<-PK_BL_topo1[!PK_BL_topo1$V1=="0",]

KK_BL_topo3<-read.delim("KatohalleriKowalyr_w100_branchlengths.topo3.txt",sep=" ",header=F)
KK_BL_topo3<-KK_BL_topo3[!KK_BL_topo3$V1=="0",]
KK_BL_topo2<-read.delim("KatohalleriKowalyr_w100_branchlengths.topo2.txt",sep=" ",header=F)
KK_BL_topo2<-KK_BL_topo2[!KK_BL_topo2$V1=="0",]
KK_BL_topo1<-read.delim("KatohalleriKowalyr_w100_branchlengths.topo1.txt",sep=" ",header=F)
KK_BL_topo1<-KK_BL_topo1[!KK_BL_topo1$V1=="0",]

BK_BL_topo3<-read.delim("BukohalleriKowalyr_w100_branchlengths.topo3.txt",sep=" ",header=F)
BK_BL_topo3<-BK_BL_topo3[!BK_BL_topo3$V1=="0",]
BK_BL_topo2<-read.delim("BukohalleriKowalyr_w100_branchlengths.topo2.txt",sep=" ",header=F)
BK_BL_topo2<-BK_BL_topo2[!BK_BL_topo2$V1=="0",]
BK_BL_topo1<-read.delim("BukohalleriKowalyr_w100_branchlengths.topo1.txt",sep=" ",header=F)
BK_BL_topo1<-BK_BL_topo1[!BK_BL_topo1$V1=="0",]

ZK_BL_topo3<-read.delim("KowahalleriZapaarenosalyr_w100_branchlengths.topo3.txt",sep=" ",header=F)
ZK_BL_topo3<-ZK_BL_topo3[!ZK_BL_topo3$V1=="0",]
ZK_BL_topo2<-read.delim("KowahalleriZapaarenosalyr_w100_branchlengths.topo2.txt",sep=" ",header=F)
ZK_BL_topo2<-ZK_BL_topo2[!ZK_BL_topo2$V1=="0",]
ZK_BL_topo1<-read.delim("KowahalleriZapaarenosalyr_w100_branchlengths.topo1.txt",sep=" ",header=F)
ZK_BL_topo1<-ZK_BL_topo1[!ZK_BL_topo1$V1=="0",]

MK<-data.frame(rbind(cbind(MK_BL_topo1$V5,"Topo1"),cbind(MK_BL_topo3$V5,"Topo3")))
colnames(MK)<-c("Branchlength","Topo")
MK$Branchlength<-as.numeric(MK$Branchlength)

shapiro.test(MK$Branchlength)
#sample size must be between 3 and 5000
hist(MK$Branchlength)
require(car)
leveneTest(MK$Branchlength~MK$Topo)
#Levene's Test for Homogeneity of Variance (center = median)
#          Df F value    Pr(>F)    
#group      1  247.25 < 2.2e-16 ***
#      116345                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Warnmeldung:
#In leveneTest.default(y = y, group = group, ...) : group coerced to factor.

t.test(MK$Branchlength[MK$Topo=="Topo1"],MK$Branchlength[MK$Topo=="Topo3"],exact=T)
#        Welch Two Sample t-test
#
#data:  MK$Branchlength[MK$Topo == "Topo1"] and MK$Branchlength[MK$Topo == "Topo3"]
#t = 100.02, df = 115847, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# 0.2128305 0.2213386
#sample estimates:
#mean of x mean of y 
#1.0777713 0.8606867

t.test(MK_BL_topo1$V5,MK_BL_topo3$V5)
#p-value < 2.2e-16
t.test(MK_BL_topo1$V5,MK_BL_topo1$V3)
#p-value < 2.2e-16
t.test(MK_BL_topo3$V5,MK_BL_topo3$V3)
#p-value < 2.2e-16


t.test(PK_BL_topo1$V5,PK_BL_topo3$V5)
#p-value < 2.2e-16
t.test(PK_BL_topo1$V5,PK_BL_topo1$V3)
#p-value < 2.2e-16
t.test(PK_BL_topo3$V5,PK_BL_topo3$V3)
#p-value < 2.2e-16


t.test(KK_BL_topo1$V5,KK_BL_topo3$V5)
#p-value < 2.2e-16
t.test(KK_BL_topo1$V5,KK_BL_topo1$V3)
#p-value < 2.2e-16
t.test(KK_BL_topo3$V5,KK_BL_topo3$V3)
#p-value < 2.2e-16


t.test(BK_BL_topo1$V5,BK_BL_topo3$V5)
#p-value < 2.2e-16
t.test(BK_BL_topo1$V5,BK_BL_topo1$V3)
#p-value < 2.2e-16
t.test(BK_BL_topo3$V5,BK_BL_topo3$V3)
#p-value < 2.2e-16

t.test(ZK_BL_topo1$V5,ZK_BL_topo3$V5)
#p-value < 2.2e-16
t.test(ZK_BL_topo1$V5,ZK_BL_topo1$V3)
#p-value < 2.2e-16
t.test(ZK_BL_topo3$V5,ZK_BL_topo3$V3)
#p-value < 2.2e-16


t.test(MK_BL_topo3$V5,BK_BL_topo3$V5)
#p-value < 2.2e-16

t.test(MK_BL_topo3$V5,KK_BL_topo3$V5)
#p-value = 0.03297

t.test(MK_BL_topo3$V5,PK_BL_topo3$V5)
#p-value < 2.2e-16

t.test(KK_BL_topo3$V5,BK_BL_topo3$V5)
#p-value < 2.2e-16

t.test(KK_BL_topo3$V5,PK_BL_topo3$V5)
#p-value < 2.2e-16

t.test(BK_BL_topo3$V5,PK_BL_topo3$V5)
#p-value < 2.2e-16







