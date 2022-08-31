data<-read.csv("str_K4_rep1_f.csv",header=FALSE,sep=";")
colnames(data)<-c("Genotype","Cluster1","Cluster2","Cluster3","Cluster4")

pdf("Structure1.pdf",width=30,height=5,pointsize=12,paper="special" )
par(mar=c(5,7,1,0),xpd=T)
par(lwd=1.3)
barplot(t(as.matrix(data[,2:5])),ylab="Assignment probability to cluster",cex.lab=1.5,cex.axis=1.3,las=2,col=rainbow(4, alpha=1),names=data$Genotype,cex.names=0.75)
dev.off()

pdf("Structure1.pdf",width=30,height=5,paper="special")
par(mar=c(5,7,1,0),xpd=T)
par(lwd=1.3)
bp<-barplot(t(as.matrix(data[-60,2:5])),ylab="Assignment probability to cluster",cex.lab=1.5,cex.axis=1.3,las=2,col=rainbow(4, alpha=1),xpd=T)
mtext(at=bp,line=0,text=as.character(data$Genotype),las=2,side=1)
dev.off()

#install.packages(c("ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)
#install.packages("rlang", type = "source")
#devtools::install_github('royfrancis/pophelper')

# load library for use
library(pophelper)
library(ggplot2)
library(gridExtra)
library(gtable)
library(tidyr)
library(grid)

qlist<-list.files(pattern="_f$")
qdata<-readQ(qlist)

library("SNPRelate")
genofile<-snpgdsOpen("F:/Introgression/Demography_redo/All5/All5_introgression_4dgsites.gds")
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

qdata_t<-tabulateQ(qdata)
qdata_s<-summariseQ(qdata_t)

qdata3<-lapply(qdata,"rownames<-",samples)

qdata_t2<-tabulateQ(as.qlist(qdata))
qdata_s2<-summariseQ(qdata_t2)
p<-evannoMethodStructure(qdata_s2,exportplot=T,returnplot=T,basesize=5,linesize=0.7,exportpath=getwd())
grid.arrange(p$plot)

slist1 <- alignK(qdata)
p1 <- plotQ(slist1,imgoutput="join",returnplot=T,exportplot=T,basesize=5,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd())
grid.arrange(p1$plot[[1]])


pdf("Strcuture1_introgression_allks.pdf",width=6,height=4,paper="special",pointsize=15)
par(mar=c(6,2,2,2))
grid.arrange(p1$plot[[1]])
dev.off()

#4 clusters
test<-qdata
grplabels<-data.frame(lab1=c(rep("A ce",3),rep("A cr",6),rep("A ly",9),rep("A pe",2),rep("A th",8),rep("Buko tetraploid",8),rep("Buko diploid",18),rep("Kato tetraploid",10),rep("Kato diploid",1),rep("Kato tetraploid",8),rep("Kowa tetraploid",8),rep("Kowa diploid",11)
,rep("Mias tetraploid",1),rep("Mias diploid",2),rep("Mias tetraploid",6),rep("Mias diploid",6),rep("Mias tetraploid",11),rep("Mias diploid",1),rep("Mias tetraploid",1),rep("Mias diploid",6),rep("Piek tetraploid",19),rep("Piek diploid",10),rep("Piek tetraploid",1)
,rep("Piek diploid",3),rep("Piek tetraploid",1),rep("Zapa diploid",12),rep("Zapa tetraploid",10),rep("Zapa diploid",4)))

slist2 <- alignK(test[28:30])

slist2 <- alignK(qdata[28:30])
p1 <- plotQ(slist2,imgoutput="join",returnplot=T,exportplot=T,basesize=3,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01)
#grid.arrange(p1$plot[[1]])

p1 <- plotQ(slist2,imgoutput="join",returnplot=T,exportplot=T,basesize=4,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=2,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9)

pdf("Structure1_introgression_k4_reordered.pdf",width=12,height=8,paper="special",pointsize=15)
par(mar=c(25,5,2,2))
par(oma=c(25,5,2,2))
par(xpd=T)
p1 <- plotQ(slist2,imgoutput="join",returnplot=T,exportplot=T,basesize=3,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01)
grid.arrange(p1$plot[[1]])
dev.off()

Mean_k4<-data.frame(samples)
Mean_k4$Cluster1<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster1,slist2$str_K4_rep2_f$Cluster1,slist2$str_K4_rep3_f$Cluster1)))
Mean_k4$Cluster2<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster2,slist2$str_K4_rep2_f$Cluster2,slist2$str_K4_rep3_f$Cluster2)))
Mean_k4$Cluster3<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster3,slist2$str_K4_rep2_f$Cluster3,slist2$str_K4_rep3_f$Cluster3)))
Mean_k4$Cluster4<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster4,slist2$str_K4_rep2_f$Cluster4,slist2$str_K4_rep3_f$Cluster4)))
colnames(Mean_k4)[1]<-"Genotype"
require(openxlsx)
write.xlsx(Mean_k4,"Structure_k4.xlsx",row.names=F,overwrite=T)


grplabels2<-data.frame(lab1=c(rep("A ce",3),rep("A cr",6),rep("A ly",9),rep("A pe",2),rep("A th",8),rep("Buko tetraploid",8),rep("Buko diploid",18),rep("Kato tetraploid",10),rep("Kato diploid",1),rep("Kato tetraploid",8),rep("Kowa tetraploid",8),rep("Kowa diploid",11)
,rep("Mias tetraploid",1),rep("Mias diploid",2),rep("Mias tetraploid",6),rep("Mias diploid",6),rep("Mias tetraploid",11),rep("Mias diploid",1),rep("Mias tetraploid",1),rep("Mias diploid",6),rep("Piek tetraploid",19),rep("Piek diploid",10),rep("Piek tetraploid",1)
,rep("Piek diploid",3),rep("Piek tetraploid",1),rep("Zapa diploid",12),rep("Zapa tetraploid",10),rep("Zapa diploid",4)))


#mean
pdf("Structure1_introgression_k4_reordered_mean.pdf",width=12,height=6,paper="special",pointsize=15)

slist3 <- mergeQ(alignK(qdata[28:30]))
p1 <- plotQ(slist3,returnplot=T,exportplot=T,basesize=4,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01,grplab=grplabels2[,1,drop=FALSE],grplabsize=5,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9)
par(mar=c(25,5,2,2))
par(oma=c(25,5,2,2))
par(xpd=T)
grid.arrange(p1$plot[[1]])
dev.off()


pdf("Structure1_introgression_k4_reordered_mean2.pdf",width=6,height=5,paper="special",pointsize=15)

slist3 <- mergeQ(alignK(qdata[28:30]))
p1 <- plotQ(slist3,returnplot=T,exportplot=T,basesize=4,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=4,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9,subset=
c("Mias diploid","Piek diploid","Buko diploid","Kato diploid","Zapa diploid","Kowa diploid","Mias tetraploid","Piek tetraploid","Buko tetraploid","Kato tetraploid","Zapa tetraploid","Kowa tetraploid","A ly","A cr","A ce","A pe","A th"))
par(mar=c(25,5,2,2))
par(oma=c(25,5,2,2))
par(xpd=T)
grid.arrange(p1$plot[[1]])

dev.off()


pdf("Structure1_introgression_k4_reordered_mean3.pdf",width=10,height=2.5,paper="special",pointsize=15)

slist3 <- mergeQ(alignK(qdata[28:30]))
p1 <- plotQ(slist3,returnplot=T,exportplot=T,basesize=4,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=4,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9,subset=
c("Mias diploid","Piek diploid","Buko diploid","Kato diploid","Zapa diploid","Kowa diploid","Mias tetraploid","Piek tetraploid","Buko tetraploid","Kato tetraploid","Zapa tetraploid","Kowa tetraploid","A ly","A cr","A ce","A pe","A th"),clustercol=c("blue","lightslateblue","red","yellow"))
par(mar=c(25,5,2,2))
par(oma=c(25,5,2,2))
par(xpd=T)
grid.arrange(p1$plot[[1]])

dev.off()

pdf("Structure1_introgression_k4_reordered_mean4.pdf",width=10,height=2.5,paper="special",pointsize=15)

slist3 <- mergeQ(alignK(qdata[28:30]))
p1 <- plotQ(slist3,returnplot=T,exportplot=T,basesize=4,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=4,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9,subset=
c("Mias diploid","Piek diploid","Buko diploid","Kato diploid","Zapa diploid","Kowa diploid","Mias tetraploid","Piek tetraploid","Buko tetraploid","Kato tetraploid","Zapa tetraploid","Kowa tetraploid","A ly","A cr","A ce","A pe","A th"),clustercol=c("darkblue","salmon2","darkgreen","bisque1"))
par(mar=c(25,5,2,2))
par(oma=c(25,5,2,2))
par(xpd=T)
grid.arrange(p1$plot[[1]])

dev.off()

Mean_k4<-data.frame(samples)
Mean_k4$Cluster1<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster1,slist2$str_K4_rep2_f$Cluster1,slist2$str_K4_rep3_f$Cluster1)))
Mean_k4$Cluster2<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster2,slist2$str_K4_rep2_f$Cluster2,slist2$str_K4_rep3_f$Cluster2)))
Mean_k4$Cluster3<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster3,slist2$str_K4_rep2_f$Cluster3,slist2$str_K4_rep3_f$Cluster3)))
Mean_k4$Cluster4<-rowMeans(as.matrix(cbind(slist2$str_K4_rep1_f$Cluster4,slist2$str_K4_rep2_f$Cluster4,slist2$str_K4_rep3_f$Cluster4)))
colnames(Mean_k4)[1]<-"Genotype"
require(openxlsx)
write.xlsx(Mean_k4,"Structure_k4.xlsx",row.names=F,overwrite=T)



#7 clusters
slist2 <- alignK(test[37:39])
p1 <- plotQ(slist2,imgoutput="join",returnplot=T,exportplot=T,basesize=4,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=2,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9)
#grid.arrange(p1$plot[[1]])

pdf("Strcuture1_introgression_k7_reordered.pdf",width=12,height=8,paper="special",pointsize=15)
par(mar=c(6,2,2,2))
grid.arrange(p1$plot[[1]])
dev.off()


#8 clusters
slist2 <- alignK(test[40:42])
p1 <- plotQ(slist2,imgoutput="join",returnplot=T,exportplot=T,basesize=4,barbordercol="black",showindlab=T,useindlab=T,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=2,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9)
#grid.arrange(p1$plot[[1]])

pdf("Strcuture1_introgression_k8_reordered.pdf",width=12,height=8,paper="special",pointsize=15)
par(mar=c(6,2,2,2))
grid.arrange(p1$plot[[1]])
dev.off()

all<-test
all2<-list()
#all, mean of 3 reps
for (i in seq(4,45,3))
	{slistk <- alignK(all[i:(i+2)])
	slistk_mean<-mergeQ(slistk)
	all2<-append(all2,slistk_mean)
	}
all3<-lapply(all2,"rownames<-",samples)


p1 <- plotQ(as.qlist(all3[c(7:14,1:6)]),imgoutput="join",returnplot=T,exportplot=F,basesize=4,barbordercol="black",showindlab=F,useindlab=F,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=2,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=300,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9,panelratio=c(1,2))
#grid.arrange(p1$plot[[1]])

pdf("Structure1_introgression_all_reordered.pdf",width=12,height=8,paper="special",pointsize=15)
par(mar=c(6,2,2,2))
grid.arrange(p1$plot[[1]])
dev.off()




shiny=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E")

p1 <- plotQ(as.qlist(all3[c(7:14,1:6)]),imgoutput="join",returnplot=T,exportplot=F,basesize=4,barbordercol="black",showindlab=F,useindlab=F,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=2,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9,clustercol=shiny[1:15],panelratio=c(1,2))

pdf("Structure1_introgression_all_reordered_diffcol.pdf",width=12,height=8,paper="special",pointsize=15)
par(mar=c(6,2,2,2))
grid.arrange(p1$plot[[1]])
dev.off()




p1a <- plotQ(as.qlist(all3[c(7:13)]),imgoutput="join",returnplot=T,exportplot=F,basesize=4,barbordercol="black",showindlab=F,useindlab=F,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=2,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9,clustercol=shiny[1:15],panelratio=c(1,1))
p1b <- plotQ(as.qlist(all3[c(13,1:6)]),imgoutput="join",returnplot=T,exportplot=F,basesize=4,barbordercol="black",showindlab=F,useindlab=F,exportpath=getwd(),panelspacer=0.01,grplab=grplabels[,1,drop=FALSE],grplabsize=2,linesize=0.8,pointsize=2,ordergrp=T,grplabangle=90,grplabheight=50,indlabangle=90,indlabvjust=1,grplabjust=1,grplabspacer=0,linepos=1,dpi=600,grplabpos=0.9,clustercol=shiny[1:15],panelratio=c(1,1))

pdf("Structure1_introgression_all_reordered_diffcol_a.pdf",width=12,height=8,paper="special",pointsize=15)
par(mar=c(6,2,2,2))
grid.arrange(p1a$plot[[1]])
dev.off()
pdf("Structure1_introgression_all_reordered_diffcol_b.pdf",width=12,height=8,paper="special",pointsize=15)
par(mar=c(6,2,2,2))
grid.arrange(p1b$plot[[1]])
dev.off()




