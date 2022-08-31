
files<-list.files(pattern="*arhacomb.dedup.MQ25.idxstats")

for (file in files)
	{nam<-gsub(".arhacomb.dedup.MQ25.idxstats","",file)
	assign(nam,read.table(file,header=F,sep="\t"))
	}

filelist<-gsub(".arhacomb.dedup.MQ25.idxstats","",files)
filelist<-filelist[c(1,9,6,7,2,5,8,10,3,4)]

for (file in filelist)
	{assign(paste("Perc_mapped_species",file,sep="_"),c())
	for (j in 1:8)
		{
		assign(paste("Perc_mapped_species",file,sep="_"),c(get(paste("Perc_mapped_species",file,sep="_")),get(file)$V3[j]/(get(file)$V3[j]+get(file)$V3[j+8])*100))
		}
	}

for (file in filelist)
	{for (j in 9:16)
		{
		assign(paste("Perc_mapped_species",file,sep="_"),c(get(paste("Perc_mapped_species",file,sep="_")),get(file)$V3[j]/(get(file)$V3[j]+get(file)$V3[j-8])*100))
		}
	}

Perc_mapped_species_list<-ls(pattern="Perc_mapped_species_*")

Perc_mapped_species_df<-matrix(unlist(lapply(Perc_mapped_species_list,function(x) get(x))),nrow=10,byrow=T)
rownames(Perc_mapped_species_df)<-Perc_mapped_species_list

require(RColorBrewer)

pdf("Heatmap_hybrids.pdf",width=15,height=20,paper="special",pointsize=20)

par(mar=c(7,10,1,1))
image(t(Perc_mapped_species_df[c(1,9,6,7,2,5,8,10,3,4),]),col=colorRampPalette(brewer.pal(9,"Greys"))(20),ylab="",xlab="",xaxt="none",yaxt="none")
axis(side=1,c(1:8,1:8),at=seq(0,1,1/15))

pops<-c(rep("Mias",1),"Zapa",rep("Mias",2),rep("Zapa",2),rep("Mias",2),rep("Zapa",2))
cols_pops<-c(rep("red",1),"black",rep("red",2),rep("black",2),rep("red",2),rep("black",2))

for (i in 1:10)
	{mtext(side=2,at=seq(0,1,length.out=10)[i],text=pops[i],las=2,line=3,font=2,adj=0,col=cols_pops[i])
	}
segments(-0.15,-0.05,-0.15,0.245, xpd = TRUE,lwd=2)
segments(-0.15,0.255,-0.15,0.645, xpd = TRUE,lwd=2)
segments(-0.15,0.655,-0.15,1.05, xpd = TRUE,lwd=2)

mtext(side=2,at=0.13,las=0,line=4,text="Potential\nHybrids",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=4,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=4,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

segments(-0.03,-0.1,0.495,-0.1, xpd = TRUE,lwd=2)
segments(0.505,-0.1,1.03,-0.1, xpd = TRUE,lwd=2)
mtext(side=1,at=0.25,las=0,line=3,text="Arabidopsis arenosa",padj=0,cex=1.3,font=3)
mtext(side=1,at=0.75,las=0,line=3,text="Arabidopsis halleri",padj=0,cex=1.3,font=3)
mtext(side=1,at=0.25,las=0,line=4,text="genome",padj=0,cex=1.3)
mtext(side=1,at=0.75,las=0,line=4,text="genome",padj=0,cex=1.3)

mtext(side=1,at=-0.2,las=0,line=1,text="Scaffold",padj=0,cex=1.3,xpd = TRUE)
mtext(side=1,at=-0.2,las=0,line=3,text="Genome",padj=0,cex=1.3,xpd = TRUE)

box()
dev.off()

pdf("Heatmap_hybrids_named.pdf",width=15,height=20,paper="special",pointsize=20)

par(mar=c(7,10,1,1))
image(t(Perc_mapped_species_df[c(1,2,10,4,5,9,11,3,6,7,8),]),col=colorRampPalette(brewer.pal(9,"Greys"))(30),ylab="",xlab="",xaxt="none",yaxt="none")
axis(side=1,c(1:8,1:8),at=seq(0,1,1/15))

pops<-c("Mias001a17","Mias001a22","Zapa002a06",rep("Mias",2),rep("Zapa",2),rep("Mias",2),rep("Zapa",2))
cols_pops<-c(rep("red",2),"black",rep("red",2),rep("black",2),rep("red",2),rep("black",2))

for (i in 1:11)
	{mtext(side=2,at=seq(0,1,length.out=11)[i],text=pops[i],las=2,line=1,font=2,adj=1,col=cols_pops[i])
	}
segments(-0.22,-0.05,-0.22,0.245, xpd = TRUE,lwd=2)
segments(-0.22,0.255,-0.22,0.645, xpd = TRUE,lwd=2)
segments(-0.22,0.655,-0.22,1.05, xpd = TRUE,lwd=2)

mtext(side=2,at=0.1,las=0,line=6.5,text="Potential\nHybrids",padj=0,cex=1.3)
mtext(side=2,at=0.48,las=0,line=6.5,text="Arabidopsis\narenosa",padj=0,cex=1.3,font=3)
mtext(side=2,at=0.85,las=0,line=6.5,text="Arabidopsis\nhalleri",padj=0,cex=1.3,font=3)

segments(-0.03,-0.1,0.495,-0.1, xpd = TRUE,lwd=2)
segments(0.505,-0.1,1.03,-0.1, xpd = TRUE,lwd=2)
mtext(side=1,at=0.25,las=0,line=3,text="Arabidopsis arenosa",padj=0,cex=1.3,font=3)
mtext(side=1,at=0.75,las=0,line=3,text="Arabidopsis halleri",padj=0,cex=1.3,font=3)
mtext(side=1,at=0.25,las=0,line=4,text="genome",padj=0,cex=1.3)
mtext(side=1,at=0.75,las=0,line=4,text="genome",padj=0,cex=1.3)

mtext(side=1,at=-0.2,las=0,line=1,text="Scaffold",padj=0,cex=1.3,xpd = TRUE)
mtext(side=1,at=-0.2,las=0,line=3,text="Genome",padj=0,cex=1.3,xpd = TRUE)

box()
dev.off()


pdf("Heatmap_hybrids_paper.pdf",width=15,height=8,paper="special",pointsize=18)
par(oma=c(1,1,1,1))
par(mar=c(7,8,1,1))
image(t(Perc_mapped_species_df[rev(c(1,9,6,7,2,5,8,10,3,4)),]),col=colorRampPalette(brewer.pal(9,"Greys"))(30),ylab="",xlab="",xaxt="none",yaxt="none")
axis(side=1,c(1:8,1:8),at=seq(0,1,1/15),cex.axis=1.3)

pops<-c("AaMias02","AaMias01","AaZapa02","AaZapa01","AhMias02","AhMias01","AhZapa02","AhZapa01","Zapa002a06","Mias001a22")
cols_pops<-rev(c(rep("red",1),"black",rep("black",2),rep("red",2),rep("black",2),rep("red",2)))

for (i in 1:10)
	{mtext(side=2,at=seq(0,1,length.out=10)[i],text=pops[i],las=2,line=1,font=2,adj=1,col=cols_pops[i],cex=1.5)
	}
segments(-0.03,-0.2,0.495,-0.2, xpd = TRUE,lwd=2.5)
segments(0.505,-0.2,1.03,-0.2, xpd = TRUE,lwd=2.5)
mtext(side=1,at=0.25,las=0,line=3,text="Arabidopsis arenosa",padj=0,cex=1.5,font=3)
mtext(side=1,at=0.75,las=0,line=3,text="Arabidopsis halleri",padj=0,cex=1.5,font=3)
mtext(side=1,at=0.25,las=0,line=4,text="genome",padj=0,cex=1.5)
mtext(side=1,at=0.75,las=0,line=4,text="genome",padj=0,cex=1.5)

mtext(side=1,at=-0.15,las=0,line=1,text="Scaffold",padj=0,cex=1.5,xpd = TRUE)
mtext(side=1,at=-0.15,las=0,line=2.5,text="Genome",padj=0,cex=1.5,xpd = TRUE)

box()
dev.off()


