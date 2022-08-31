
CNV_up<-read.table("CNVs_fisher.txt",sep="\t",header=T)
GS<-read.table("GS_fisher.txt",sep="\t",header=T)
Twisst<-read.table("Twisst_fisher.txt",sep="\t",header=T)



bpCNV_fisher_sig<-barplot(sort(-log10(CNV$q_value)),las=3,xlim=c(0,max(-log10(CNV$q_value))),width=1.2)
bpGS_fisher_sig<-barplot(sort(-log10(GS$q_value)),las=3,xlim=c(0,max(-log10(GS$q_value))),width=1.2)
bpTwisst_fisher_sig<-barplot(sort(-log10(Twisst$q_value)),las=3,xlim=c(0,max(-log10(Twisst$q_value))),width=1.2)



pdf('GOenrichment_Mapman_CNV_GS_Twisst_qvalue.pdf', width=9, height=9,paper="special",pointsize=15)
par(oma=c(6,8,2,1))
layout(1:3,heights=c(1,1,0.75))
par(mar=c(2,15,2,0),xpd=T)

barplot(sort(-log10(CNV$q_value)),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,8),col=rgb(0.1,0.4,1,alpha=0.3),cex.axis=1.3,axes=F,main="Copy number divergence",cex.main=1.5)
axis(side=2,tick=F,labels=rev(c("Stress","Mitochondrial electron transport/\n ATP synthesis","Nutrition&ionome","Metal\nhomeostasis")),at=bpCNV_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)

barplot(sort(-log10(GS$q_value)),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,8),col=rgb(1,0.2,0.2,alpha=0.3),cex.axis=1.3,axes=F,main="Selection",cex.main=1.5)
axis(side=2,tick=F,labels=rev(c("Nutrition&ionome","Metal\nhomeostasis","Transport","DNA")),at=bpGS_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(-log[10](adjusted~italic(p)-value)),side=1,cex=1.3,outer=T,line=1)

barplot(sort(-log10(Twisst$q_value)),beside=T,horiz=T,las=1,cex.lab=1.5,xlim=c(0,8),width=1.2,col=rgb(0.2,1,0.2,alpha=0.3),cex.axis=1.5,main="Introgression",cex.main=1.5)
axis(side=2,tick=F,labels=rev(c("Nutrition&ionome","Transport","Metal\nhomeostasis")),at=bpTwisst_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(Enriched~Mapman~categories),side=2,cex=1.3,outer=T,line=4,xlim=c(0,8))

dev.off()


bpCNV_fisher_sig<-barplot(sort(CNV$Fold_enrichment),las=3,xlim=c(0,max(CNV$Fold_enrichment)),width=1.2)
bpGS_fisher_sig<-barplot(sort(GS$Fold_enrichment),las=3,xlim=c(0,max(GS$Fold_enrichment)),width=1.2)
bpTwisst_fisher_sig<-barplot(sort(Twisst$Fold_enrichment),las=3,xlim=c(0,max(Twisst$Fold_enrichment)),width=1.2)



pdf('GOenrichment_Mapman_CNV_GS_Twisst_Fold_enrichment.pdf', width=9, height=9,paper="special",pointsize=15)
par(oma=c(6,8,2,1))
layout(1:3,heights=c(1,1,0.75))
par(mar=c(2,15,2,0),xpd=T)

barplot(sort(CNV$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,8),col=rgb(0.1,0.4,1,alpha=0.3),cex.axis=1.3,axes=F,main="Copy number divergence",cex.main=1.5)
axis(side=2,tick=F,labels=rev(c("Mitochondrial electron transport/\n ATP synthesis","Nutrition&ionome","Metal\nhomeostasis","Stress")),at=bpCNV_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)

barplot(sort(GS$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,8),col=rgb(1,0.2,0.2,alpha=0.3),cex.axis=1.3,axes=F,main="Selection",cex.main=1.5)
axis(side=2,tick=F,labels=rev(c("Nutrition&ionome","Metal\nhomeostasis","DNA","Transport")),at=bpGS_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(-log[10](adjusted~italic(p)-value)),side=1,cex=1.3,outer=T,line=1)

barplot(sort(Twisst$Fold_enrichment),beside=T,horiz=T,las=1,cex.lab=1.5,xlim=c(0,8),width=1.2,col=rgb(0.2,1,0.2,alpha=0.3),cex.axis=1.5,main="Introgression",cex.main=1.5)
axis(side=2,tick=F,labels=rev(c("Nutrition&ionome","Metal\nhomeostasis","Transport")),at=bpTwisst_fisher_sig,xpd=T,srt=45,las=2,font=2,cex.axis=1.3,cex=1.3)
mtext(text=expression(Enriched~Mapman~categories),side=2,cex=1.3,outer=T,line=4,xlim=c(0,8))

dev.off()




