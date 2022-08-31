options(java.parameters = "-Xmx12g")
require(openxlsx)

####################################
#Folded AF plots with Krom and Kosi#
####################################
M_K<-read.delim("MiasarenosaK.table",header=T,sep="\t",quote="")
M_K<-M_K[,1:4]
K_M<-read.delim("KowaarenosaM.table",header=T,sep="\t",quote="")
K_M<-K_M[,1:4]

M_ha<-read.delim("Miasarenosaha.table",header=T,sep="\t",quote="")
M_ha<-M_ha[,1:4]
ha_M<-read.delim("haM.table",header=T,sep="\t",quote="")
ha_M<-ha_M[,1:4]

AFdata_Aa<-M_K[,1:2]
AFdata_Aa$AF_M<-M_K$AC/M_K$AN
AFdata_Aa$AF_K<-K_M$AC/K_M$AN
names(AFdata_Aa)<-c("Chrom","POS","AF_M","AF_K")

AFdata_Mh<-M_ha[,1:2]
AFdata_Mh$AF_M<-M_ha$AC/M_ha$AN
AFdata_Mh$AF_ha<-ha_M$AC/ha_M$AN
names(AFdata_Mh)<-c("Chrom","POS","AF_M","AF_ha")



repol<-read.table("F:/Introgression/Demography_redo/All5/Introgression_07repolarized.lookupKey.minInd_4.txt",sep="\t",header=F)
remove<-read.table("F:/Introgression/Demography_redo/All5/Repol_exclude.txt",sep="\t",header=F)
AFdata_Aa<-AFdata_Aa[!(paste0(AFdata_Aa$Chrom,AFdata_Aa$POS)%in%paste0(remove$V1,remove$V2)),]
repol2<-repol[paste0(repol$V1,repol$V2)%in%paste0(AFdata_Aa$Chrom,AFdata_Aa$POS),]
#373483 out of of 734892
AFdata_Aa2<-AFdata_Aa
AFdata_Aa2$AF_M<-ifelse(paste0(AFdata_Aa2$Chrom,AFdata_Aa2$POS)%in%paste0(repol2$V1,repol2$V2),1-AFdata_Aa2$AF_M,AFdata_Aa2$AF_M)
AFdata_Aa2$AF_K<-ifelse(paste0(AFdata_Aa2$Chrom,AFdata_Aa2$POS)%in%paste0(repol2$V1,repol2$V2),1-AFdata_Aa2$AF_K,AFdata_Aa2$AF_K)
AFdata_Aa<-AFdata_Aa2


AFdata_Mh<-AFdata_Mh[!(paste0(AFdata_Mh$Chrom,AFdata_Mh$POS)%in%paste0(remove$V1,remove$V2)),]
repol2<-repol[paste0(repol$V1,repol$V2)%in%paste0(AFdata_Mh$Chrom,AFdata_Mh$POS),]
#355685 out of of 734892
AFdata_Mh2<-AFdata_Mh
AFdata_Mh2$AF_M<-ifelse(paste0(AFdata_Mh2$Chrom,AFdata_Mh2$POS)%in%paste0(repol2$V1,repol2$V2),1-AFdata_Mh2$AF_M,AFdata_Mh2$AF_M)
AFdata_Mh2$AF_ha<-ifelse(paste0(AFdata_Mh2$Chrom,AFdata_Mh2$POS)%in%paste0(repol2$V1,repol2$V2),1-AFdata_Mh2$AF_ha,AFdata_Mh2$AF_ha)
AFdata_Mh<-AFdata_Mh2

#################
#Candidate genes#
#################
Candidates2<-read.xlsx("Selected_genes_Fst_introgression_1_MiasKowa.xlsx",1)

Candidates5<-Candidates2[!duplicated(Candidates2$Alyr_ID),]
Candidates3<-Candidates5[!is.na(Candidates5$Alyr_ID),]
Candidates6<-Candidates3[!Candidates3$Alyr_ID=="",]

genelist<-read.table("F:/Cu_project/Lyrata/Alyrata_384_v2.1.gene.gff3",sep="\t",header=F)
colnames(genelist)<-c("Scaffold","Source","Type","gene_start","gene_end","Score","Strand","Phase","Attributes")
genelist<-genelist[genelist$Type=="gene",]
ID<-substr(as.character(genelist$Attributes),4,nchar(as.character(genelist$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
genelist$Alyr_ID<-ID3[,1]

Candidates<-merge(Candidates6,genelist,by="Alyr_ID")
Candidates<-Candidates[!is.na(Candidates$Name),]

genemiddle=list((genelist$gene_start[match(Candidates$Alyr_ID,genelist$Alyr_ID)]+genelist$gene_end[match(Candidates$Alyr_ID,genelist$Alyr_ID)])/2)
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
#windowsize<-((Candidates$gene_end-Candidates$gene_start)+1)*2 #bp genemiddle +/- windowsize will be displayed

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest_Aa", j, sep ="")
 	AF1_Aa<-AFdata_Aa[AFdata_Aa$Chrom==Candidates$Chr_gene[j],]
	AFofinterest_Aa<-AF1_Aa[AF1_Aa$POS>=(genemiddle[[1]][j]-windowsize)&AF1_Aa$POS<=(genemiddle[[1]][j]+windowsize),]
	assign(nam, AFofinterest_Aa)
	}
rm(AFofinterest_Aa)
list<-ls(pattern="^AFofinterest_Aa")

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest_Mh", j, sep ="")
 	AF1_Mh<-AFdata_Mh[AFdata_Mh$Chrom==Candidates$Chr_gene[j],]
	AFofinterest_Mh<-AF1_Mh[AF1_Mh$POS>=(genemiddle[[1]][j]-windowsize)&AF1_Mh$POS<=(genemiddle[[1]][j]+windowsize),]
	assign(nam, AFofinterest_Mh)
	}
rm(AFofinterest_Mh)
list<-ls(pattern="^AFofinterest_Mh")

for (j in 1:nrow(Candidates))
	{nam <- paste("Ath_IDsofinterest_Aa", j, sep ="")
 	Ath_IDs1_Aa<-genelist[genelist$Scaffold==Candidates$Chr_gene[j],]
	Ath_IDsofinterest_Aa<-Ath_IDs1_Aa[Ath_IDs1_Aa$gene_start>=(genemiddle[[1]][j]-windowsize)&Ath_IDs1_Aa$gene_start<=(genemiddle[[1]][j]+windowsize)|Ath_IDs1_Aa$gene_end>=(genemiddle[[1]][j]-windowsize)&Ath_IDs1_Aa$gene_end<=(genemiddle[[1]][j]+windowsize),]
	assign(nam, Ath_IDsofinterest_Aa)
	}
rm(Ath_IDsofinterest_Aa)

test<-as.factor(AFdata_Aa$Chrom)
AFdata_Aa[,1]<-test

winsize=10
AFmatrix_Aa=data.frame()
for (i in levels(AFdata_Aa[,1]))
        {datanow<-AFdata_Aa[AFdata_Aa[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,4])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix_Aa=rbind(AFmatrix_Aa,wmatrix)
        }
names(AFmatrix_Aa)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Mias","AF_Kowa")


test<-as.character(AFdata_Aa$Chrom)
AFdata_Aa[,1]<-test

test2<-AFmatrix_Aa
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix_Aa[,i]))
	}
AFmatrix_Aa<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest_Aa", j, sep ="")
 	AF1_Aa<-AFmatrix_Aa[AFmatrix_Aa$Chrom==Candidates$Chr_gene[j],]
	AFfreqofinterest_Aa<-AF1_Aa[AF1_Aa$Window_start>=(genemiddle[[1]][j]-windowsize-100000)&AF1_Aa$Window_start<=(genemiddle[[1]][j]+windowsize+100000)|AF1_Aa$Window_end>=(genemiddle[[1]][j]-windowsize-100000)&AF1_Aa$Window_end<=(genemiddle[[1]][j]+windowsize+100000),]
	assign(nam, AFfreqofinterest_Aa)
	}
rm(AFfreqofinterest_Aa)




test<-as.factor(AFdata_Mh$Chrom)
AFdata_Mh[,1]<-test

winsize=10
AFmatrix_Mh=data.frame()
for (i in levels(AFdata_Mh[,1]))
        {datanow<-AFdata_Mh[AFdata_Mh[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,4])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix_Mh=rbind(AFmatrix_Mh,wmatrix)
        }
names(AFmatrix_Mh)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Mias","AF_halleri")


test<-as.character(AFdata_Mh$Chrom)
AFdata_Mh[,1]<-test

test2<-AFmatrix_Mh
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix_Mh[,i]))
	}
AFmatrix_Mh<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest_Mh", j, sep ="")
 	AF1_Mh<-AFmatrix_Mh[AFmatrix_Mh$Chrom==Candidates$Chr_gene[j],]
	AFfreqofinterest_Mh<-AF1_Mh[AF1_Mh$Window_start>=(genemiddle[[1]][j]-windowsize-100000)&AF1_Mh$Window_start<=(genemiddle[[1]][j]+windowsize+100000)|AF1_Mh$Window_end>=(genemiddle[[1]][j]-windowsize-100000)&AF1_Mh$Window_end<=(genemiddle[[1]][j]+windowsize+100000),]
	assign(nam, AFfreqofinterest_Mh)
	}
rm(AFfreqofinterest_Mh)

#SNPeffdf<-read.table("GS_KromKosi/SNPeff_KromKosi.table",sep="\t",header=T)

#AFdatasnp<-cbind(AFdata[,1:3],AFdata[,4])
#names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
#SNPeffdfAF<-merge(SNPeffdf,AFdatasnp)
#SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
#SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
#SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]

genelist2<-read.table("F:/Cu_project/Lyrata/Alyrata_384_v2.1.gene_exons.gff3",sep="\t",header=F)
colnames(genelist2)<-c("Scaffold","Source","Type","gene_start","gene_end","Score","Strand","Phase","Attributes")
exonlist<-na.exclude(genelist2[genelist2$Type=="exon",])
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
ID<-substr(as.character(exonlist$Attributes),4,nchar(as.character(exonlist$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
exonlist$Alyr_ID<-ID3[,1]
exonlist<-exonlist[!duplicated(exonlist[,-9]),]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-na.exclude(exonlist[exonlist$Scaffold==Candidates$Chr_gene[j],])
	Exonsofinterest<-Exons1[Exons1$gene_start>=(genemiddle[[1]][j]-windowsize)&Exons1$gene_start<=(genemiddle[[1]][j]+windowsize)|Exons1$gene_end>=(genemiddle[[1]][j]-windowsize)&Exons1$gene_end<=(genemiddle[[1]][j]+windowsize),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-na.exclude(exonlist[exonlist$Scaffold==Candidates$Chr_gene[j],])
	Exonsofinterestcand<-Exons1[Exons1$Alyr_ID==Candidates$Alyr_ID[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)

Candidates<-merge(Candidates,genelist,by="Alyr_ID")
arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$Strand.x=="+",2,1)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Ath_IDsofinterest_Aa",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Denser_AF_plot_MiasKowa_HMgenes_",j,"_",Candidates$Alyr_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=8, height=6,paper="special")
	par(mar=c(5,5,1,1)+0.1)
	par(oma=c(0,0,1,10)+0.1)
	par(mgp=c(3,0.7,0))
	plot((get(paste("AFofinterest_Aa",j,sep=""))$AF_M-get(paste("AFofinterest_Aa",j,sep=""))$AF_K)~get(paste("AFofinterest_Aa",j,sep=""))$POS,ylab="",ylim=c(-1,1.2),cex.lab=1,xlab="",xlim=c(genemiddle[[1]][j]-windowsize,genemiddle[[1]][j]+windowsize),xaxt='n',yaxt='n',col="grey30")
	rect(Candidates$gene_start.x[j],-2, Candidates$gene_end.x[j],2,col="lightgrey",border = NA)	
	points((get(paste("AFofinterest_Aa",j,sep=""))$AF_M-get(paste("AFofinterest_Aa",j,sep=""))$AF_K)~get(paste("AFofinterest_Aa",j,sep=""))$POS,col="grey30")
	text(x=genemiddle[[1]][j],y=1.15,labels=Candidates$Name[j],col="red")
	for(i in 1:nrow(get(paste("Ath_IDsofinterest_Aa",j,sep="")))) 
		{arrows(x0=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=1.1,x1=Candidates$gene_end.x[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	arrowdir3<-ifelse(get(paste("Exonsofinterestcand",j,sep=""))$Strand=="+",2,1)
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir3[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest_Aa",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)],lwd=2,col="black")
	lines(get(paste("AFfreqofinterest_Aa",j,sep=""))$AF_Mias[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)],lwd=2,col="red")
	lines(get(paste("AFfreqofinterest_Mh",j,sep=""))$AF_halleri[order(get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle)],lwd=2,col="blue")

	legend("topright",legend=c(expression(25~SNP-window~bold(Mias)),expression(25~SNP-window~bold(Kowa)),expression(25~SNP-window~bold(halleri))),pch=15,col=c("red","black","blue"),bg="white",bty="n",xpd=NA,inset=c(-0.59,0.02),title="Average allele frequency")
#	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])],col="purple")
#	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])],col="orange")
#	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])],col="red")
#	legend("topright",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2,title="Predicted effect",bty="n",xpd=NA,inset=c(-0.35,0.2),title.adj=2)
	axis(side=1,cex=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex.axis=0.8,adj=0.1)
	axis(side=2,cex=1,at=axTicks(side=2),labels=axTicks(side=2),las=2,cex.axis=0.8)
	box()
	mtext(text="Scaffold position (kb)",side=1,line=2,cex=1)
	mtext(text=expression(Allele~frequency~difference~bold(Mias)~-~bold(Kowa)),side=2,line=2,cex=1)
	dev.off()
	}




for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Ath_IDsofinterest_Aa",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Denser_Diff_AF_plot_MiasKowahalleri_HMgenes_",j,"_",Candidates$Alyr_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=7, height=7,paper="special")
	par(mfrow=c(2,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(3,0,1,1)+0.1)
	par(mgp=c(3,0.7,0))
	plot((get(paste("AFofinterest_Aa",j,sep=""))$AF_M-get(paste("AFofinterest_Aa",j,sep=""))$AF_K)~get(paste("AFofinterest_Aa",j,sep=""))$POS,ylab="",ylim=c(-1,1.2),cex.lab=1,xlab="",xlim=c(genemiddle[[1]][j]-windowsize,genemiddle[[1]][j]+windowsize),xaxt='n',yaxt='n',col="grey30")
	rect(Candidates$gene_start.x[j],-2, Candidates$gene_end.x[j],2,col="lightgrey",border = NA)	
	points((get(paste("AFofinterest_Aa",j,sep=""))$AF_M-get(paste("AFofinterest_Aa",j,sep=""))$AF_K)~get(paste("AFofinterest_Aa",j,sep=""))$POS,col="grey30")
	text(x=genemiddle[[1]][j],y=1.15,labels=Candidates$Name[j],col="red")
	for(i in 1:nrow(get(paste("Ath_IDsofinterest_Aa",j,sep="")))) 
		{arrows(x0=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=1.1,x1=Candidates$gene_end.x[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	arrowdir3<-ifelse(get(paste("Exonsofinterestcand",j,sep=""))$Strand=="+",2,1)
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir3[i],length=0,col="red",lwd=3)
		}
#	lines(get(paste("AFfreqofinterest_Aa",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)]~
#	get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)],lwd=2,col="black")
#	lines(get(paste("AFfreqofinterest_Aa",j,sep=""))$AF_Mias[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)]~
#	get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest_Aa",j,sep=""))$Window_middle)],lwd=2,col="red")
#	legend("topright",legend=c(expression(25~SNP-window~bold(Mias)),expression(25~SNP-window~bold(Kowa)),expression(25~SNP-window~bold(halleri))),pch=15,col=c("red","black","blue"),bg="white",bty="n",xpd=NA,inset=c(-0.59,0.02),title="Average allele frequency")
#	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])],col="purple")
#	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])],col="orange")
#	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])],col="red")
#	legend("topright",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2,title="Predicted effect",bty="n",xpd=NA,inset=c(-0.35,0.2),title.adj=2)
	axis(side=1,cex=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex.axis=0.8,adj=0.1)
	axis(side=2,cex=1,at=axTicks(side=2),labels=axTicks(side=2),las=2,cex.axis=0.8)
	mtext(text=expression(Allele~frequency~bold(Mias)~-~bold(Kowa)),side=2,line=2,cex=1)
	plot((get(paste("AFofinterest_Mh",j,sep=""))$AF_M-get(paste("AFofinterest_Mh",j,sep=""))$AF_ha)~get(paste("AFofinterest_Mh",j,sep=""))$POS,ylab="",ylim=c(-1,1.2),cex.lab=1,xlab="",xlim=c(genemiddle[[1]][j]-windowsize,genemiddle[[1]][j]+windowsize),xaxt='n',yaxt='n',col="grey30")
	rect(Candidates$gene_start.x[j],-2, Candidates$gene_end.x[j],2,col="lightgrey",border = NA)	
	points((get(paste("AFofinterest_Mh",j,sep=""))$AF_M-get(paste("AFofinterest_Mh",j,sep=""))$AF_ha)~get(paste("AFofinterest_Mh",j,sep=""))$POS,col="grey30")
	text(x=genemiddle[[1]][j],y=1.15,labels=Candidates$Name[j],col="red")
	for(i in 1:nrow(get(paste("Ath_IDsofinterest_Aa",j,sep="")))) 
		{arrows(x0=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=1.1,x1=Candidates$gene_end.x[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	arrowdir3<-ifelse(get(paste("Exonsofinterestcand",j,sep=""))$Strand=="+",2,1)
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir3[i],length=0,col="red",lwd=3)
		}
#	lines(get(paste("AFfreqofinterest_Mh",j,sep=""))$AF_Mias[order(get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle)]~
#	get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle)],lwd=2,col="red")
#	lines(get(paste("AFfreqofinterest_Mh",j,sep=""))$AF_halleri[order(get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle)]~
#	get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest_Mh",j,sep=""))$Window_middle)],lwd=2,col="blue")

#	legend("topright",legend=c(expression(25~SNP-window~bold(Mias)),expression(25~SNP-window~bold(Kowa)),expression(25~SNP-window~bold(halleri))),pch=15,col=c("red","black","blue"),bg="white",bty="n",xpd=NA,inset=c(-0.7,0.02),title="Average allele frequency")
#	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Alyr_ID[j])],col="purple")
#	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Alyr_ID[j])],col="orange")
#	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Alyr_ID[j])],col="red")
#	legend("topright",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2,title="Predicted effect",bty="n",xpd=NA,inset=c(-0.35,0.2),title.adj=2)
	axis(side=1,cex=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex.axis=0.8,adj=0.1)
	axis(side=2,cex=1,at=axTicks(side=2),labels=axTicks(side=2),las=2,cex.axis=0.8)
	box()
	mtext(text="Scaffold position (kb)",side=1,line=1,cex=1,outer=T)
	mtext(text=expression(Allele~frequency~bold(Mias)~-~bold(halleri)),side=2,line=2,cex=1)


	dev.off()
	}



for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Ath_IDsofinterest_Aa",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Denser_Diff_AF_plot_MiasKowahalleri_HMgenes_",j,"_",Candidates$Alyr_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=13, height=9,paper="special",pointsize=15)
	par(mfrow=c(2,1))
	par(mar=c(1,5,0,1)+0.1)
	par(oma=c(3,1,1,1)+0.1)
	par(mgp=c(3,0.7,0))
	plot((get(paste("AFofinterest_Aa",j,sep=""))$AF_M-get(paste("AFofinterest_Aa",j,sep=""))$AF_K)~get(paste("AFofinterest_Aa",j,sep=""))$POS,ylab="",ylim=c(-1,1.2),cex.lab=1,xlab="",xlim=c(genemiddle[[1]][j]-windowsize,genemiddle[[1]][j]+windowsize),xaxt='n',yaxt='n',col="black",pch=19)
	rect(genemiddle[[1]][j]-windowsize-100000,-0.5, genemiddle[[1]][j]+windowsize+100000,0.5,col=rgb(1,0.2,0.2,alpha=0.3),border = NA)	
	rect(Candidates$gene_start.x[j],-2, Candidates$gene_end.x[j],2,col=rgb(0,0,0,alpha=0.3),border = NA)	
	points((get(paste("AFofinterest_Aa",j,sep=""))$AF_M-get(paste("AFofinterest_Aa",j,sep=""))$AF_K)~get(paste("AFofinterest_Aa",j,sep=""))$POS,col="black",pch=19)
#	text(x=genemiddle[[1]][j],y=1.15,labels=Candidates$Name[j],col="red")
	for(i in 1:nrow(get(paste("Ath_IDsofinterest_Aa",j,sep="")))) 
		{arrows(x0=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=2.5)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=1.1,x1=Candidates$gene_end.x[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=2.5)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
		}
	arrowdir3<-ifelse(get(paste("Exonsofinterestcand",j,sep=""))$Strand=="+",2,1)
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir3[i],length=0,col="red",lwd=5)
		}
	axis(side=1,cex=1,at=axTicks(side=1),labels=NA,cex.axis=1.3,adj=0.1)
	axis(side=2,cex=1,at=axTicks(side=2),labels=axTicks(side=2),las=2,cex.axis=1.3)
#	mtext(text=expression(Allele~frequency~bold(Mias)~-~bold(Kowa)),side=2,line=3,cex=1.5,cex.lab=1.5)
	plot((get(paste("AFofinterest_Mh",j,sep=""))$AF_M-get(paste("AFofinterest_Mh",j,sep=""))$AF_ha)~get(paste("AFofinterest_Mh",j,sep=""))$POS,ylab="",ylim=c(-1,1.2),cex.lab=1,xlab="",xlim=c(genemiddle[[1]][j]-windowsize,genemiddle[[1]][j]+windowsize),xaxt='n',yaxt='n',col="black",pch=19)
	rect(genemiddle[[1]][j]-windowsize-100000,-0.5, genemiddle[[1]][j]+windowsize+100000,0.5,col=rgb(1,0.2,0.2,alpha=0.3),border = NA)	
	rect(Candidates$gene_start.x[j],-2, Candidates$gene_end.x[j],2,col=rgb(0,0,0,alpha=0.3),border = NA)	
	points((get(paste("AFofinterest_Mh",j,sep=""))$AF_M-get(paste("AFofinterest_Mh",j,sep=""))$AF_ha)~get(paste("AFofinterest_Mh",j,sep=""))$POS,col="black",pch=19)
#	text(x=genemiddle[[1]][j],y=1.15,labels=Candidates$Name[j],col="red")
	for(i in 1:nrow(get(paste("Ath_IDsofinterest_Aa",j,sep="")))) 
		{arrows(x0=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Ath_IDsofinterest_Aa",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey40",lwd=2.5)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=1.1,x1=Candidates$gene_end.x[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=2.5)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=1.1,code=1,length=0,col="grey40",lwd=5)
		}
	arrowdir3<-ifelse(get(paste("Exonsofinterestcand",j,sep=""))$Strand=="+",2,1)
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=1.1,code=arrowdir3[i],length=0,col="red",lwd=5)
		}
	axis(side=1,cex=1,at=axTicks(side=1),labels=NA,cex.axis=1.3,adj=0.1)
	axis(side=2,cex=1,at=axTicks(side=2),labels=axTicks(side=2),las=2,cex.axis=1.3)
	box()
#	mtext(text="Scaffold position (Mb)",side=1,line=1,cex=1.5,outer=T,cex.lab=1.5)
#	mtext(text=expression(Allele~frequency~bold(Mias)~-~bold(halleri)),side=2,line=3,cex=1.5,cex.lab=1.5)


	dev.off()
	}

#formatC(axTicks(side=1)/1000000,big.mark=",",format="f",digits=2)
