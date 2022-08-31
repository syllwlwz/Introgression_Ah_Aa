require(openxlsx)
Species<-read.xlsx("Halleri_arenosa_flow_cytometry_formerge.xlsx",1)

Flow<-read.csv("Flow_cytometry.csv",sep=";")

All<-merge(Species,Flow,by="Sample_ID",all.y=T)

write.xlsx(All,"Sequ_flow_cytometry_complete.xlsx")


pdf("Zn_Field_plant_flow_cyt.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,5,1,1))
boxplot(log10(All$Zn)~All$Ploidy*All$Population,las=2,ylim=c(0,5),col=c("darkgreen","darkblue","darkcyan","purple"),names=c(rep("",15)),xlab="",ylab=expression(Log[10](Leaf~Zn*" [ppm]")),xaxt="n",cex.lab=1.5,cex.axis=1.3)
stripchart(log10(All$Zn)~All$Ploidy*All$Population,vertical=TRUE,add=TRUE, pch=21,bg="black")
abline(h=log10(3000),col="red")
segments(0.6,-0.25, 4.4,-0.25, xpd = TRUE)
segments(4.6,-0.25, 8.4,-0.25, xpd = TRUE)
segments(8.6,-0.25, 12.4,-0.25, xpd = TRUE)
segments(12.6,-0.25, 15.4,-0.25, xpd = TRUE)
mtext(c("Mias","Piek","Kato","Buko"),side=1,line=2,at=c(2.5,6.5,10.5,14),cex=1.5,col=c("red"))
par(lwd=2)
legend("topright",fill=c("darkgreen","darkblue","darkcyan","purple"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa),italic(Intermediate),italic(Lyrata~intermediate))),cex=1.5,bty="n")
dev.off()










list<-read.table("D:/Introgression/ICP_arenosa2/Sequencing_list_introgression_final.csv",sep=",",header=T)


data<-read.xlsx("D:/Introgression/ICP_arenosa2/Field_halleri_arenosa_withVero_dedup.xlsx",1)
Flow2<-read.csv("Flow_cytometry_Species.csv",sep=";")
Flow1<-rbind(Flow,Flow2)
Flow1<-Flow1[!(duplicated(Flow1)&(Flow1$Sample_ID!="Piek_h6"&Flow1$Sample_ID!="Piek_h10")),]

icp<-merge(data,Flow1,by="Sample_ID",all.y=T)
#icp<-icp[!icp$Population=="Kowa",]

write.xlsx(icp,"Sequ_flow_cytometry_Vero.xlsx")

icp$Population<-factor(icp$Population,levels=c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"))

pdf("Zn_Field_plant_flow_cyt_SPP.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,5,1,1))
boxplot(log10(icp$Zn)~icp$Ploidy*icp$Population,las=2,ylim=c(0,5),col=c("darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Zn*" [ppm]")),xaxt="n",cex.lab=1.5,cex.axis=1.3)
stripchart(log10(icp$Zn)~icp$Ploidy*icp$Population,vertical=TRUE,add=TRUE, pch=21,bg="black")
abline(h=log10(3000),col="red")
mtext(c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"),side=1,line=2,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),cex=1.5,col=c("red","black","red","black","red","red","red"))
abline(v=2.5,col="grey")
abline(v=4.5,col="grey")
abline(v=6.5,col="grey")
abline(v=8.5,col="grey")
abline(v=10.5,col="grey")
abline(v=12.5,col="grey")
par(lwd=2)
legend("bottomleft",fill=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bg="white",box.lwd=0.5)
dev.off()

pdf("Cd_Field_plant_flow_cyt.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,5,1,1))
boxplot(log10(icp$Cd)~icp$Ploidy*icp$Population,las=2,ylim=c(-1,4),border=c("darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Cd)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
abline(h=log10(100),col="red")
mtext(c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"),side=1,line=2,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),cex=1.5,col=c("red","black","red","black","red","red","red"))
abline(v=2.5,col="grey")
abline(v=4.5,col="grey")
abline(v=6.5,col="grey")
abline(v=8.5,col="grey")
abline(v=10.5,col="grey")
abline(v=12.5,col="grey")
par(lwd=2)
legend("topleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bg="white",box.lwd=0.5)
dev.off()

pdf("Pb_Field_plant_flow_cyt.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,5,1,1))
boxplot(log10(icp$Pb)~icp$Ploidy*icp$Population,las=2,ylim=c(-2,4),border=c("darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Pb)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
abline(h=log10(1000),col="red")
mtext(c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"),side=1,line=2,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),cex=1.5,col=c("red","black","red","black","red","red","red"))
abline(v=2.5,col="grey")
abline(v=4.5,col="grey")
abline(v=6.5,col="grey")
abline(v=8.5,col="grey")
abline(v=10.5,col="grey")
abline(v=12.5,col="grey")
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bg="white",box.lwd=0.5)
dev.off()

soil<-read.xlsx("ICP/Soil_icp_wVero_wRicardo.xlsx")
All<-merge(icp,soil,by="Sample_ID")
All$Population<-factor(All$Population.x,levels=c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"))

pdf("Zn_Field_plant_flow_cyt_ratio.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,5,1,1))
Zn<-log10(All$Zn/All$CH_Zn)
boxplot(Zn~All$Ploidy*All$Population.x,las=2,ylim=c(0,5),col=c("darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Zn/Exchangeable~Zn*" [ppm]")),xaxt="n",cex.lab=1.5,cex.axis=1.3)
stripchart(Zn~All$Ploidy*All$Population.x,vertical=TRUE,add=TRUE, pch=21,bg="black")
mtext(c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"),side=1,line=2,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),cex=1.5,col=c("red","black","red","black","red","red","red"))
abline(v=2.5,col="grey")
abline(v=4.5,col="grey")
abline(v=6.5,col="grey")
abline(v=8.5,col="grey")
abline(v=10.5,col="grey")
abline(v=12.5,col="grey")
par(lwd=2)
legend("topright",fill=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bg="white",box.lwd=0.5)
dev.off()

pdf("Cd_Field_plant_flow_cyt_ratio.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,5,1,1))
Cd<-log10(All$Cd/(All$CH_Cd+1))
boxplot(Cd~All$Ploidy*All$Population.x,las=2,ylim=c(-1,3),border=c("darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Cd/(Exchangeable~Cd+1))*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
mtext(c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"),side=1,line=2,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),cex=1.5,col=c("red","black","red","black","red","red","red"))
abline(v=2.5,col="grey")
abline(v=4.5,col="grey")
abline(v=6.5,col="grey")
abline(v=8.5,col="grey")
abline(v=10.5,col="grey")
abline(v=12.5,col="grey")
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bg="white",box.lwd=0.5)
dev.off()

pdf("Pb_Field_plant_flow_cyt_ratio.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,5,1,1))
Pb<-log10(All$Pb/(All$CH_Pb+1))
boxplot(Pb~All$Ploidy*All$Population.x,las=2,ylim=c(0,3),border=c("darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Pb/(Exchangeable~Pb+1))*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
mtext(c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"),side=1,line=2,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),cex=1.5,col=c("red","black","red","black","red","red","red"))
abline(v=2.5,col="grey")
abline(v=4.5,col="grey")
abline(v=6.5,col="grey")
abline(v=8.5,col="grey")
abline(v=10.5,col="grey")
abline(v=12.5,col="grey")
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bg="white",box.lwd=0.5)
dev.off()

#load Ricardo M NM halleri for comp
Ricardo<-read.xlsx("Leaf and soil data_import.xlsx",2)
Zn_R<-log10(Ricardo$Leaf_Zn/Ricardo$Exchangeable_Zn)
Zn_M<-Zn_R[Ricardo$Soil_cluster=="Metalliferous"]
Zn_NM<-Zn_R[Ricardo$Soil_cluster=="Non_metalliferous"]

pdf("Zn_Field_plant_flow_cyt_ratio_wRicardo.pdf",width=8,height=8,paper="special",pointsize=14)
layout(mat = matrix(c(1,2),nrow = 1,ncol = 2),heights = c(2, 2),widths = c(3, 1))
par(lwd=1)
par(mar=c(6,5,1,1))
Zn<-log10(All$Zn/All$CH_Zn)
x<-boxplot(Zn~All$Ploidy*All$Population.x,las=2,ylim=c(0,6),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Zn/Exchangeable~Zn)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
y<-boxplot(list(Zn_M,Zn_NM),las=2,ylim=c(0,6),names=c(rep("",2)),xlab="",ylab=expression(Log[10](Leaf~Zn/Exchangeable~Zn)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3,yaxt="n")
boxplot(x$stats[,c(5:10,13:14,1:4,11:12)],las=2,ylim=c(0,6),border=c("darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Zn)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
mtext(c("Klet","Kowa","Mias","Zapa","Buko","Kato","Piek"),side=1,line=2,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),cex=1.5,col=c("red","black","red","black","red","red","red"))
abline(v=2.5,col="grey")
abline(v=4.5,col="grey")
abline(v=6.5,col="grey")
abline(v=8.5,col="grey")
abline(v=10.5,col="grey")
abline(v=12.5,col="grey")
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bg="white",box.lwd=0.5)

boxplot(y$stats,las=2,ylim=c(0,6),border=c("darkgreen"),names=c(rep("",2)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,yaxt="n")



dev.off()




