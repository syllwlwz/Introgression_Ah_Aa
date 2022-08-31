
require(openxlsx)
###############
#Metrics plots#
###############
Candidates2<-read.xlsx("Selected_genes_Fst_introgression_1_MiasKowa.xlsx",1)

Candidates5<-Candidates2[!duplicated(Candidates2$Alyr_ID),]
Candidates3<-Candidates5[!is.na(Candidates5$Alyr_ID),]
Candidates6<-Candidates3[!Candidates3$Alyr_ID=="",]

options(java.parameters = "-Xmx12g")

filelist = list.files(pattern = glob2rx("Fst_Mias_Kowa_arenosa_scaffold*_25SNPwins.table")) 
require(gtools)
filelist <- mixedsort(filelist)
fst_25 = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t"))) #merge all count files

#################################
#kill windows bigger than 100 kb#
#################################
fst_25<-fst_25[fst_25$End_position-fst_25$Start_position<100000,]

genelist<-read.table("F:/Cu_project/Lyrata/Alyrata_384_v2.1.gene.gff3",sep="\t",header=F)
colnames(genelist)<-c("Scaffold","Source","Type","gene_start","gene_end","Score","Strand","Phase","Attributes")
genelist<-genelist[genelist$Type=="gene",]
ID<-substr(as.character(genelist$Attributes),4,nchar(as.character(genelist$Attributes)))
ID2<-lapply(strsplit(ID,"\\."),"[",1)
ID3<-data.frame(matrix(unlist(ID2),nrow=length(ID2),byrow=T))
genelist$Alyr_ID<-ID3[,1]

Candidates<-merge(Candidates6,genelist,by="Alyr_ID")
Candidates<-Candidates[!is.na(Candidates$Name),]

Candidates$gene_size<-abs(Candidates$gene_end.x-Candidates$gene_start.x)
genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-50000 #bp genemiddle +/- windowsize will be displayed
test3<-as.character(Candidates$Scaffold)
Candidates$Chr<-test3

AFdataIntro<-fst_25

winmiddleIntro<-c((AFdataIntro$Start_pos+AFdataIntro$End_pos)/2)
AFdataIntro<-cbind(AFdataIntro,winmiddleIntro)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterestIntro", j, sep ="")
 	AF1Intro<-AFdataIntro[tolower(AFdataIntro$Scaffold)==Candidates$Chr[j],]
	AFofinterestIntro<-AF1Intro[AF1Intro$winmiddleIntro>=(Candidates$gene_start.x[j]-windowsize-100000)&AF1Intro$winmiddleIntro<=(Candidates$gene_end.x[j]+windowsize+100000),]
	assign(nam, AFofinterestIntro)
	}
rm(AFofinterestIntro)
list<-ls(pattern="^AFofinterestIntro")

genemiddle=list((Candidates$gene_start.x+Candidates$gene_end.x)/2)

#Quantiles:
Fst_percIntro<-quantile(AFdataIntro$Fst_win,0.995)

for (j in 1:nrow(Candidates))
	{nam <- paste("Alyr_IDsofinterest", j, sep ="")
 	Alyr_IDs1<-genelist[tolower(genelist$Scaffold)==Candidates$Chr[j],]
	Alyr_IDsofinterest<-Alyr_IDs1[Alyr_IDs1$gene_start>=(genemiddle[[1]][j]-windowsize)&Alyr_IDs1$gene_start<=(genemiddle[[1]][j]+windowsize)|Alyr_IDs1$gene_end>=(genemiddle[[1]][j]-windowsize)&Alyr_IDs1$gene_end<=(genemiddle[[1]][j]+windowsize),]
	assign(nam, Alyr_IDsofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$Strand=="+",2,1)

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
 	Exons1<-na.exclude(exonlist[tolower(exonlist$Scaffold)==Candidates$Chr[j],])
	Exonsofinterest<-Exons1[Exons1$gene_start>=(genemiddle[[1]][j]-windowsize)&Exons1$gene_start<=(genemiddle[[1]][j]+windowsize)|Exons1$gene_end>=(genemiddle[[1]][j]-windowsize)&Exons1$gene_end<=(genemiddle[[1]][j]+windowsize),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$Strand=="+",2,1)

Candidates$Alyr_ID<-as.character(Candidates$Alyr_ID)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$Alyr_ID==Candidates$Alyr_ID[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)

d<-Candidates


filelist2 = list.files(pattern = glob2rx("NeisD_Mias_arenosa_allhalleri_scaffold*_25SNPwins.table")) 
require(gtools)
filelist2 <- mixedsort(filelist2)
NeisDhalleri_25 = do.call(rbind,lapply(filelist2,function(fn)read.delim(fn,header=T, sep="\t"))) #merge all count files
NeisDhalleri_25<-NeisDhalleri_25[(NeisDhalleri_25$End_position-NeisDhalleri_25$Start_position)<100000,]

AFdataIntrohalleri<-NeisDhalleri_25
winmiddleIntro<-c((AFdataIntrohalleri$Start_pos+AFdataIntrohalleri$End_pos)/2)
AFdataIntrohalleri<-cbind(AFdataIntrohalleri,winmiddleIntro)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterestIntrohalleri", j, sep ="")
 	AF1Introhalleri<-AFdataIntrohalleri[tolower(AFdataIntrohalleri$Scaffold)==Candidates$Chr[j],]
	AFofinterestIntrohalleri<-AF1Introhalleri[AF1Introhalleri$winmiddleIntro>=(Candidates$gene_start.x[j]-windowsize[j]-100000)&AF1Introhalleri$winmiddleIntro<=(Candidates$gene_end.x[j]+windowsize[j]+100000),]
	assign(nam, AFofinterestIntrohalleri)
	}
rm(AFofinterestIntrohalleri)
list<-ls(pattern="^AFofinterestIntrohalleri")
NeisD_percIntrohalleri<-quantile(AFdataIntrohalleri$NeisD,0.005)

Dxy_MK<-read.delim("DxyMiasKowaarenosa.csv",sep="\t",header=F)
colnames(Dxy_MK)<-c("Scaffold","Start_pos","End_pos","Dxy")
Dxy_MK_25<-Dxy_MK[(Dxy_MK$End_pos-Dxy_MK$Start_pos)<100000,]

AFdataDxy_MK<-Dxy_MK_25
winmiddleDxy_MK<-c((AFdataDxy_MK$Start_pos+AFdataDxy_MK$End_pos)/2)
AFdataDxy_MK<-cbind(AFdataDxy_MK,winmiddleDxy_MK)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterestDxy_MK", j, sep ="")
 	AF1Dxy_MK<-AFdataDxy_MK[tolower(AFdataDxy_MK$Scaffold)==Candidates$Chr[j],]
	AFofinterestDxy_MK<-AF1Dxy_MK[AF1Dxy_MK$winmiddleDxy_MK>=(Candidates$gene_start.x[j]-windowsize[j]-100000)&AF1Dxy_MK$winmiddleDxy_MK<=(Candidates$gene_end.x[j]+windowsize[j]+100000),]
	assign(nam, AFofinterestDxy_MK)
	}
rm(AFofinterestDxy_MK)
list<-ls(pattern="^AFofinterestDxy_MK")
Dxy_percMK<-quantile(AFdataDxy_MK$Dxy,0.995,na.rm=TRUE)

Dxy_Mha<-read.delim("DxyMiashalleri.csv",sep="\t",header=F)
colnames(Dxy_Mha)<-c("Scaffold","Start_pos","End_pos","Dxy")
Dxy_Mha_25<-Dxy_Mha[(Dxy_Mha$End_pos-Dxy_Mha$Start_pos)<100000,]

AFdataDxy_Mha<-Dxy_Mha_25
winmiddleDxy_Mha<-c((AFdataDxy_Mha$Start_pos+AFdataDxy_Mha$End_pos)/2)
AFdataDxy_Mha<-cbind(AFdataDxy_Mha,winmiddleDxy_Mha)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterestDxy_Mha", j, sep ="")
 	AF1Dxy_Mha<-AFdataDxy_Mha[tolower(AFdataDxy_Mha$Scaffold)==Candidates$Chr[j],]
	AFofinterestDxy_Mha<-AF1Dxy_Mha[AF1Dxy_Mha$winmiddleDxy_Mha>=(Candidates$gene_start.x[j]-windowsize[j]-100000)&AF1Dxy_Mha$winmiddleDxy_Mha<=(Candidates$gene_end.x[j]+windowsize[j]+100000),]
	assign(nam, AFofinterestDxy_Mha)
	}
rm(AFofinterestDxy_Mha)
list<-ls(pattern="^AFofinterestDxy_Mha")
Dxy_percMha<-quantile(AFdataDxy_Mha$Dxy,0.005,na.rm=TRUE)

for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Alyr_IDsofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Mias_Fst_Dxy_Kowa_halleri",j,"_",Candidates$Alyr_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=8, height=8,paper="special")
	par(mfcol=c(2,1))
	par(mar=c(0,2,1,1)+0.1)
	par(oma=c(5,5,3,0))
	par(mgp=c(5,1,0))


	plot(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,ylab="",ylim=c(min(AFdataIntro$Fst_win),max(AFdataIntro$Fst_win)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start.x[j],min(AFdataIntro$Fst_win)-0.7,Candidates$gene_end.x[j],max(AFdataIntro$Fst_win)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,lwd=2,col="red")
	abline(h=Fst_percIntro,lty=2,col="black")
	text(x=genemiddle[[1]][j],y=1.4,labels=Candidates$Name[j],col="red",font=3)
	for(i in 1:nrow(get(paste("Alyr_IDsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=max(AFdataIntro$Fst_win)+0.2,x1=Candidates$gene_end.x[j],y1=max(AFdataIntro$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestDxy_Mha",j,sep=""))$Dxy~get(paste("AFofinterestDxy_Mha",j,sep=""))$winmiddleDxy_Mha,ylab="",ylim=c(min(AFdataDxy_Mha$Dxy,na.rm=T),max(AFdataDxy_Mha$Dxy,na.rm=T)+0.3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start.x[j],min(AFdataDxy_Mha$Dxy,na.rm=T)-0.7,Candidates$gene_end.x[j],max(AFdataDxy_Mha$Dxy,na.rm=T)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestDxy_Mha",j,sep=""))$Dxy~get(paste("AFofinterestDxy_Mha",j,sep=""))$winmiddleDxy_Mha,lwd=2,col="blue")
	abline(h=Dxy_percMha,lty=2,col="black")
	text(x=genemiddle[[1]][j],y=0.73,labels=Candidates$Name[j],col="red",font=3)
	for(i in 1:nrow(get(paste("Alyr_IDsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,x1=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,x1=Candidates$gene_end.x[j],y1=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=max(AFdataDxy_Mha$Dxy,na.rm=T)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)
	mtext(text=expression(bold(Scaffold~position~(kb))),side=1,line=4,outer=TRUE,cex=1.5)

	mtext(text=expression(bold(F[ST])~Mias-Kowa~italic(A.a.)),side=2,line=1,outer=TRUE,cex=1.5,adj=0.9,col="red")
	mtext(text=expression(bold(Dxy)~Mias~italic(A.a.-A.h.)),side=2,line=1,outer=TRUE,cex=1.5,adj=0.1,col="blue")

	mtext(text=expression(bold("Mias")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.5)
	
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)

	dev.off()
	}



for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Alyr_IDsofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Mias_Fst_NeisD_Kowa_halleri",j,"_",Candidates$Alyr_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=8, height=8,paper="special")
	par(mfcol=c(2,1))
	par(mar=c(0,2,1,1)+0.1)
	par(oma=c(5,5,3,0))
	par(mgp=c(5,1,0))


	plot(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,ylab="",ylim=c(min(AFdataIntro$Fst_win),max(AFdataIntro$Fst_win)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start.x[j],min(AFdataIntro$Fst_win)-0.7,Candidates$gene_end.x[j],max(AFdataIntro$Fst_win)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,lwd=2,col="red")
	abline(h=Fst_percIntro,lty=2,col="black")
	text(x=genemiddle[[1]][j],y=1.4,labels=Candidates$Name[j],col="red",font=3)
	for(i in 1:nrow(get(paste("Alyr_IDsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=max(AFdataIntro$Fst_win)+0.2,x1=Candidates$gene_end.x[j],y1=max(AFdataIntro$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestIntrohalleri",j,sep=""))$NeisD~get(paste("AFofinterestIntrohalleri",j,sep=""))$winmiddleIntro,ylab="",ylim=c(min(AFdataIntrohalleri$NeisD),0.8),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start.x[j],min(AFdataIntrohalleri$NeisD)-0.7,Candidates$gene_end.x[j],max(AFdataIntrohalleri$NeisD)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestIntrohalleri",j,sep=""))$NeisD~get(paste("AFofinterestIntrohalleri",j,sep=""))$winmiddleIntro,lwd=2,col="blue")
	abline(h=NeisD_percIntrohalleri,lty=2,col="black")
	text(x=genemiddle[[1]][j],y=0.73,labels=Candidates$Name[j],col="red",font=3)
	for(i in 1:nrow(get(paste("Alyr_IDsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_start[i],y0=0.7,x1=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_end[i],y1=0.7,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=0.7,x1=Candidates$gene_end.x[j],y1=0.7,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=0.7,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=0.7,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=0.7,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=0.7,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)
	mtext(text=expression(bold(Scaffold~position~(kb))),side=1,line=4,outer=TRUE,cex=1.5)

	mtext(text=expression(bold(F[ST])~Mias-Kowa~italic(A.a.)),side=2,line=1,outer=TRUE,cex=1.5,adj=0.9,col="red")
	mtext(text=expression(bold(Neis~D)~Mias~italic(A.a.-A.h.)),side=2,line=1,outer=TRUE,cex=1.5,adj=0.1,col="blue")

	mtext(text=expression(bold("Mias")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.5)
	
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)

	dev.off()
	}



















	plot(get(paste("AFofinterestDxy_MK",j,sep=""))$Dxy~get(paste("AFofinterestDxy_MK",j,sep=""))$winmiddleDxy_MK,ylab="",ylim=c(min(AFdataDxy_MK$Dxy,na.rm=T),max(AFdataDxy_MK$Dxy,na.rm=T)+0.3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start.x[j],min(AFdataDxy_MK$Dxy,na.rm=T)-0.7,Candidates$gene_end.x[j],max(AFdataDxy_MK$Dxy,na.rm=T)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestDxy_MK",j,sep=""))$Dxy~get(paste("AFofinterestDxy_MK",j,sep=""))$winmiddleDxy_MK,lwd=2,col="blue")
	abline(h=Dxy_percMK,lty=2,col="black")
	text(x=genemiddle[[1]][j],y=0.73,labels=Candidates$Name[j],col="red",font=3)
	for(i in 1:nrow(get(paste("Alyr_IDsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,x1=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,x1=Candidates$gene_end.x[j],y1=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=max(AFdataDxy_MK$Dxy,na.rm=T)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Alyr_IDsofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Mias_Fst_only_Kowa_halleri",j,"_",Candidates$Alyr_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=8, height=4,paper="special")
	par(mar=c(0,2,1,1)+0.1)
	par(oma=c(5,5,3,0))
	par(mgp=c(5,1,0))


	plot(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,ylab="",ylim=c(min(AFdataIntro$Fst_win),max(AFdataIntro$Fst_win)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start.x[j],min(AFdataIntro$Fst_win)-0.7,Candidates$gene_end.x[j],max(AFdataIntro$Fst_win)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,lwd=2,col="red")
	abline(h=Fst_percIntro,lty=2,col="black")
	text(x=genemiddle[[1]][j],y=1.4,labels=Candidates$Name[j],col="red",font=3)
	for(i in 1:nrow(get(paste("Alyr_IDsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start.x[j],y0=max(AFdataIntro$Fst_win)+0.2,x1=Candidates$gene_end.x[j],y1=max(AFdataIntro$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)
	mtext(text=expression(bold(Scaffold~position~(kb))),side=1,line=4,outer=TRUE,cex=1.5)

	mtext(text=expression(bold(F[ST])~Mias-Kowa~italic(A.a.)),side=2,line=1,outer=TRUE,cex=1.5,adj=0.5,col="red")
	mtext(text=expression(bold("Mias")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.5)
	
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)

	dev.off()
	}


j=10
for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Alyr_IDsofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Mias_Fst_paper_Kowa_halleri",j,"_",Candidates$Alyr_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=15, height=9,paper="special",pointsize=22)
	par(mar=c(6,6,1,1)+0.1)
	par(oma=c(1,1,1,0))
	par(mgp=c(3,0.75,0))


	plot(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,ylab="",ylim=c(0,1),cex=1.5,cex.lab=1.5,cex.axis=1.3,type="l",lwd=3,col="red",xlim=c(Candidates$gene_start.x[j]-windowsize,Candidates$gene_end.x[j]+windowsize),xaxt="n",xlab="",las=1)
	rect(Candidates$gene_start.x[j],min(AFdataIntro$Fst_win)-0.7,Candidates$gene_end.x[j],max(AFdataIntro$Fst_win)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestIntro",j,sep=""))$Fst_win~get(paste("AFofinterestIntro",j,sep=""))$winmiddleIntro,lwd=3,col="red")
	abline(h=Fst_percIntro,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=1.4,labels=Candidates$Name[j],col="red",font=3)
#	for(i in 1:nrow(get(paste("Alyr_IDsofinterest",j,sep="")))) 
#		{arrows(x0=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Alyr_IDsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
#		}
#	arrows(x0=Candidates$gene_start.x[j],y0=max(AFdataIntro$Fst_win)+0.2,x1=Candidates$gene_end.x[j],y1=max(AFdataIntro$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
#	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
#		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=1,length=0,col="grey",lwd=3)
#		}
#	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
#		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$gene_start[i],y0=max(AFdataIntro$Fst_win)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$gene_end[i],y1=max(AFdataIntro$Fst_win)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
#		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000000,big.mark=",",format="f",preserve.width="common",digits=2),cex=1.2,cex.lab=1.2,cex.axis=1.3)
	mtext(text=expression(Scaffold~position~(Mb)),side=1,line=3,cex=1.5)

	mtext(text=expression(F[ST]),side=2,line=3,cex=1.5,adj=0.5,col="black")
	box()
	
	dev.off()
	}


