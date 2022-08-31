#Spider Charts 
#replace 0 with next min/2 because of division by min

#including Vero's ans Ricardo'sdata for Mias and Zapa
require(openxlsx)
icp<-read.xlsx("nQuire_evaluation_df_updated_wKato.xlsx",1)

xtabs(~icp$Species+icp$Population)

#replace 0 with min/2
Mins<-c()
for (i in c(37:54,61:76,81:97,102:120))
	{icp[,i]<-as.numeric(icp[,i])
	Mins[i]<-min(icp[,i][icp[,i]!=0])
	}
MIns<-Mins/2

for (i in c(37:54,61:76,81:97,102:120))
	{for (j in 1:nrow(icp))
		{icp[j,i]<-ifelse(icp[j,i]==0,MIns[i],icp[j,i])
		}
	}

require(fmsb)

#adjust levels
icp$Population <- factor(icp$Population)
icp$Species <- factor(icp$Species)

#16 elements  <- c(Pb, Cd, Zn, Mn, Ni,  K, S,   P, B,  Cr, Al,  Mg, Ca, Fe, Cu, H+) the graph is made COUNTERCLOCK-wise!!

CH_col <- c( which( colnames(icp)=="CH_Pb" ), which( colnames(icp)=="CH_Cd" ), which( colnames(icp)=="CH_Zn" ),which( colnames(icp)=="CH_Mn" ),
             which( colnames(icp)=="CH_K" ),which( colnames(icp)=="CH_S" ),which( colnames(icp)=="CH_P" ),which( colnames(icp)=="CH_B" ),which( colnames(icp)=="CH_Cr" ),
            which( colnames(icp)=="CH_Al" ),which( colnames(icp)=="CH_Ca" ),which( colnames(icp)=="CH_Fe" ),which( colnames(icp)=="CH_Cu" ),which( colnames(icp)=="H" ))

names(icp[CH_col])

Max<-log10(max(as.numeric(apply(icp[,CH_col],2,function(x) {x/min(x,na.rm=T)})),na.rm=T))


###################################################################################
#Buko arenosa
infoACH <- icp[icp$Population=="Buko"&icp$Species=="arenosa",1:2]
ArawCH <- icp[icp$Population=="Buko"&icp$Species=="arenosa",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_BukoCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_BukoCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Buko",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_BukoCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Buko_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Buko_Aa <- data.frame(t(Buko_Aa))
names(Buko_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Buko halleri
infoACH <- icp[icp$Population=="Buko"&icp$Species=="halleri",1:2]
ArawCH <- icp[icp$Population=="Buko"&icp$Species=="halleri",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_BukoCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_BukoCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Buko",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_BukoCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Buko_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Buko_Ah <- data.frame(t(Buko_Ah))
names(Buko_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")


#Kato arenosa
infoACH <- icp[icp$Population=="Kato"&icp$Species=="arenosa",1:2]
ArawCH <- icp[icp$Population=="Kato"&icp$Species=="arenosa",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KatoCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_KatoCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Kato",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KatoCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Kato_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kato_Aa <- data.frame(t(Kato_Aa))
names(Kato_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Kato halleri
infoACH <- icp[icp$Population=="Kato"&icp$Species=="halleri",1:2]
ArawCH <- icp[icp$Population=="Kato"&icp$Species=="halleri",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KatoCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_KatoCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Kato",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KatoCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Kato_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kato_Ah <- data.frame(t(Kato_Ah))
names(Kato_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")



#Mias arenosa
infoACH <- icp[icp$Population=="Mias"&icp$Species=="arenosa",1:2]
ArawCH <- icp[icp$Population=="Mias"&icp$Species=="arenosa",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_MiasCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_MiasCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Mias",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_MiasCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Mias_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Mias_Aa <- data.frame(t(Mias_Aa))
names(Mias_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Mias halleri
infoACH <- icp[icp$Population=="Mias"&icp$Species=="halleri",1:2]
ArawCH <- icp[icp$Population=="Mias"&icp$Species=="halleri",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_MiasCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_MiasCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Mias",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_MiasCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Mias_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Mias_Ah <- data.frame(t(Mias_Ah))
names(Mias_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")


#Piek arenosa
infoACH <- icp[icp$Population=="Piek"&icp$Species=="arenosa",1:2]
ArawCH <- icp[icp$Population=="Piek"&icp$Species=="arenosa",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_PiekCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_PiekCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Piek",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_PiekCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Piek_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Piek_Aa <- data.frame(t(Piek_Aa))
names(Piek_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Piek halleri
infoACH <- icp[icp$Population=="Piek"&icp$Species=="halleri",1:2]
ArawCH <- icp[icp$Population=="Piek"&icp$Species=="halleri",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_PiekCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_PiekCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Piek",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_PiekCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Piek_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Piek_Ah <- data.frame(t(Piek_Ah))
names(Piek_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")



#Kowa arenosa
infoACH <- icp[icp$Population=="Kowa"&icp$Species=="arenosa",1:2]
ArawCH <- icp[icp$Population=="Kowa"&icp$Species=="arenosa",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KowaCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_KowaCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Kowa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KowaCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Kowa_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kowa_Aa <- data.frame(t(Kowa_Aa))
names(Kowa_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Kowa halleri
infoACH <- icp[icp$Population=="Kowa"&icp$Species=="halleri",1:2]
ArawCH <- icp[icp$Population=="Kowa"&icp$Species=="halleri",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KowaCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_KowaCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Kowa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KowaCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Kowa_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kowa_Ah <- data.frame(t(Kowa_Ah))
names(Kowa_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")


#Zapa arenosa
infoACH <- icp[icp$Population=="Zapa"&icp$Species=="arenosa",1:2]
ArawCH <- icp[icp$Population=="Zapa"&icp$Species=="arenosa",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_ZapaCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_ZapaCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Zapa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_ZapaCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)


Zapa_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Zapa_Aa <- data.frame(t(Zapa_Aa))
names(Zapa_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Zapa halleri
infoACH <- icp[icp$Population=="Zapa"&icp$Species=="halleri",1:2]
ArawCH <- icp[icp$Population=="Zapa"&icp$Species=="halleri",CH_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_ZapaCH       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawCH)){
  
  stats_ZapaCH  [,i-7] <- quantile(ArawCH[infoACH$Population=="Zapa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_ZapaCH  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawCH)
head(ArawCH)

As <- data.frame(Pb = ((ArawCH[,1]/min(icp[,CH_col][,1], na.rm=T))),   Cd = ((ArawCH[,2]/min(icp[,CH_col][,2], na.rm=T))),     Zn = ((ArawCH[,3]/min(icp[,CH_col][,3], na.rm=T))),
                 Mn = ((ArawCH[,4]/min(icp[,CH_col][,4], na.rm=T))),   K = ((ArawCH[,5]/min(icp[,CH_col][,5], na.rm=T))),
                 S = ((ArawCH[,6]/min(icp[,CH_col][,6], na.rm=T))),    P = ((ArawCH[,7]/min(icp[,CH_col][,7], na.rm=T))),      B = ((ArawCH[,8]/min(icp[,CH_col][,8], na.rm=T))),
                 Cr = ((ArawCH[,9]/min(icp[,CH_col][,9], na.rm=T))), Al = ((ArawCH[,10]/min(icp[,CH_col][,10], na.rm=T))),   Ca = ((ArawCH[,11]/min(icp[,CH_col][,11], na.rm=T))),
			Fe = ((ArawCH[,12]/min(icp[,CH_col][,12], na.rm=T))),   Cu = ((ArawCH[,13]/min(icp[,CH_col][,13], na.rm=T))),
                 H = ((ArawCH[,14]/min(icp[,CH_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoACH,A)

Zapa_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Zapa_Ah <- data.frame(t(Zapa_Ah))
names(Zapa_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

pdf("CH_Introgression_minall_species_separated_wKato.pdf",width=20,height=8,paper="special",pointsize=14)
layout(matrix(1:12, ncol = 6))
par(mar = c(0, 0.5,1,0.5))
par(oma = c(0,4,2,0))

dat<- radarchart(Kowa_Aa, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkblue",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Kowary),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Kowa_Ah, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkgreen",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Zapa_Aa, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkblue",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Zakopane),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Zapa_Ah, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkgreen",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Kato_Aa, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkblue",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Katowice),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Kato_Ah, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkgreen",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Buko_Aa, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkblue",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Bukowno),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Buko_Ah, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkgreen",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Piek_Aa, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkblue",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Piekary~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Piek_Ah, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkgreen",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Mias_Aa, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkblue",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Miastezko~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Mias_Ah, axistype = 1, axislabcol="grey40", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("darkgreen",5)), pty=32, plty=c(3,2,1,2,3))

dev.off()


pdf("CH_Introgression_minall_species_separated_wKato2.pdf",width=20,height=8,paper="special",pointsize=14)
layout(matrix(1:12, ncol = 6))
par(mar = c(0, 0.5,1,0.5))
par(oma = c(0,4,2,0))

dat<- radarchart(rbind(Kowa_Aa[1:2,],Kowa_Ah[5,],Kowa_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Kowary),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Kowa_Ah[1:2,],Kowa_Aa[5,],Kowa_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Zapa_Aa[1:2,],Zapa_Ah[5,],Zapa_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Zakopane),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Zapa_Ah[1:2,],Zapa_Aa[5,],Zapa_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Kato_Aa[1:2,],Kato_Ah[5,],Kato_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Katowice),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Kato_Ah[1:2,],Kato_Aa[5,],Kato_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Buko_Aa[1:2,],Buko_Ah[5,],Buko_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Bukowno),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Buko_Ah[1:2,],Buko_Aa[5,],Buko_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Piek_Aa[1:2,],Piek_Ah[5,],Piek_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Piekary~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Piek_Ah[1:2,],Piek_Aa[5,],Piek_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Mias_Aa[1:2,],Mias_Ah[5,],Mias_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Miastezko~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Mias_Ah[1:2,],Mias_Aa[5,],Mias_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dev.off()



#################################
#Extractable soil concentrations#
#################################

X_col <- c( which( colnames(icp)=="X_Pb" ), which( colnames(icp)=="X_Cd" ), which( colnames(icp)=="X_Zn" ),which( colnames(icp)=="X_Mn" ),
             which( colnames(icp)=="X_K" ),which( colnames(icp)=="X_S" ),which( colnames(icp)=="X_P" ),which( colnames(icp)=="X_B" ),which( colnames(icp)=="X_Cr" ),
            which( colnames(icp)=="X_Al" ),which( colnames(icp)=="X_Ca" ),which( colnames(icp)=="X_Fe" ),which( colnames(icp)=="X_Cu" ),which( colnames(icp)=="H" ))

names(icp[X_col])
Max<-log10(max(as.numeric(apply(icp[,X_col],2,function(x) {x/min(x,na.rm=T)})),na.rm=T))

###################################################################################
#Buko arenosa
infoAX <- icp[icp$Population=="Buko"&icp$Species=="arenosa",1:2]
ArawX <- icp[icp$Population=="Buko"&icp$Species=="arenosa",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_BukoX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_BukoX  [,i-7] <- quantile(ArawX[infoAX$Population=="Buko",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_BukoX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Buko_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Buko_Aa <- data.frame(t(Buko_Aa))
names(Buko_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Buko halleri
infoAX <- icp[icp$Population=="Buko"&icp$Species=="halleri",1:2]
ArawX <- icp[icp$Population=="Buko"&icp$Species=="halleri",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_BukoX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_BukoX  [,i-7] <- quantile(ArawX[infoAX$Population=="Buko",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_BukoX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Buko_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Buko_Ah <- data.frame(t(Buko_Ah))
names(Buko_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")


#Kato arenosa
infoAX <- icp[icp$Population=="Kato"&icp$Species=="arenosa",1:2]
ArawX <- icp[icp$Population=="Kato"&icp$Species=="arenosa",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KatoX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_KatoX  [,i-7] <- quantile(ArawX[infoAX$Population=="Kato",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KatoX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Kato_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kato_Aa <- data.frame(t(Kato_Aa))
names(Kato_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Kato halleri
infoAX <- icp[icp$Population=="Kato"&icp$Species=="halleri",1:2]
ArawX <- icp[icp$Population=="Kato"&icp$Species=="halleri",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KatoX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_KatoX  [,i-7] <- quantile(ArawX[infoAX$Population=="Kato",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KatoX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Kato_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kato_Ah <- data.frame(t(Kato_Ah))
names(Kato_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")



#Mias arenosa
infoAX <- icp[icp$Population=="Mias"&icp$Species=="arenosa",1:2]
ArawX <- icp[icp$Population=="Mias"&icp$Species=="arenosa",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_MiasX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_MiasX  [,i-7] <- quantile(ArawX[infoAX$Population=="Mias",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_MiasX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Mias_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Mias_Aa <- data.frame(t(Mias_Aa))
names(Mias_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Mias halleri
infoAX <- icp[icp$Population=="Mias"&icp$Species=="halleri",1:2]
ArawX <- icp[icp$Population=="Mias"&icp$Species=="halleri",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_MiasX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_MiasX  [,i-7] <- quantile(ArawX[infoAX$Population=="Mias",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_MiasX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Mias_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Mias_Ah <- data.frame(t(Mias_Ah))
names(Mias_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")


#Piek arenosa
infoAX <- icp[icp$Population=="Piek"&icp$Species=="arenosa",1:2]
ArawX <- icp[icp$Population=="Piek"&icp$Species=="arenosa",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_PiekX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_PiekX  [,i-7] <- quantile(ArawX[infoAX$Population=="Piek",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_PiekX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Piek_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Piek_Aa <- data.frame(t(Piek_Aa))
names(Piek_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Piek halleri
infoAX <- icp[icp$Population=="Piek"&icp$Species=="halleri",1:2]
ArawX <- icp[icp$Population=="Piek"&icp$Species=="halleri",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_PiekX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_PiekX  [,i-7] <- quantile(ArawX[infoAX$Population=="Piek",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_PiekX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Piek_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Piek_Ah <- data.frame(t(Piek_Ah))
names(Piek_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")



#Kowa arenosa
infoAX <- icp[icp$Population=="Kowa"&icp$Species=="arenosa",1:2]
ArawX <- icp[icp$Population=="Kowa"&icp$Species=="arenosa",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KowaX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_KowaX  [,i-7] <- quantile(ArawX[infoAX$Population=="Kowa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KowaX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Kowa_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kowa_Aa <- data.frame(t(Kowa_Aa))
names(Kowa_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Kowa halleri
infoAX <- icp[icp$Population=="Kowa"&icp$Species=="halleri",1:2]
ArawX <- icp[icp$Population=="Kowa"&icp$Species=="halleri",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_KowaX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_KowaX  [,i-7] <- quantile(ArawX[infoAX$Population=="Kowa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KowaX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Kowa_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kowa_Ah <- data.frame(t(Kowa_Ah))
names(Kowa_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")


#Zapa arenosa
infoAX <- icp[icp$Population=="Zapa"&icp$Species=="arenosa",1:2]
ArawX <- icp[icp$Population=="Zapa"&icp$Species=="arenosa",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_ZapaX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_ZapaX  [,i-7] <- quantile(ArawX[infoAX$Population=="Zapa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_ZapaX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Zapa_Aa <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Zapa_Aa <- data.frame(t(Zapa_Aa))
names(Zapa_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")

#Zapa halleri
infoAX <- icp[icp$Population=="Zapa"&icp$Species=="halleri",1:2]
ArawX <- icp[icp$Population=="Zapa"&icp$Species=="halleri",X_col]    

#X quantiles
#this part is only for saving quantiles statistics
stats_ZapaX       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawX)){
  
  stats_ZapaX  [,i-7] <- quantile(ArawX[infoAX$Population=="Zapa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_ZapaX  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawX)
head(ArawX)

As <- data.frame(Pb = ((ArawX[,1]/min(icp[,X_col][,1], na.rm=T))),   Cd = ((ArawX[,2]/min(icp[,X_col][,2], na.rm=T))),     Zn = ((ArawX[,3]/min(icp[,X_col][,3], na.rm=T))),
                 Mn = ((ArawX[,4]/min(icp[,X_col][,4], na.rm=T))),   K = ((ArawX[,5]/min(icp[,X_col][,5], na.rm=T))),
                 S = ((ArawX[,6]/min(icp[,X_col][,6], na.rm=T))),    P = ((ArawX[,7]/min(icp[,X_col][,7], na.rm=T))),      B = ((ArawX[,8]/min(icp[,X_col][,8], na.rm=T))),
                 Cr = ((ArawX[,9]/min(icp[,X_col][,9], na.rm=T))), Al = ((ArawX[,10]/min(icp[,X_col][,10], na.rm=T))),   Ca = ((ArawX[,11]/min(icp[,X_col][,11], na.rm=T))),
			Fe = ((ArawX[,12]/min(icp[,X_col][,12], na.rm=T))),   Cu = ((ArawX[,13]/min(icp[,X_col][,13], na.rm=T))),
                 H = ((ArawX[,14]/min(icp[,X_col][,14], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoAX,A)

Zapa_Ah <- data.frame(max = rep(Max,14), min= c(rep(0,14)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Zapa_Ah <- data.frame(t(Zapa_Ah))
names(Zapa_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu","H+")


pdf("X_Introgression_minall_species_separated_wKato.pdf",width=20,height=8,paper="special")

layout(matrix(1:12, ncol = 6))
layout(matrix(1:12, ncol = 6))
par(mar = c(0, 0.5,1,0.5))
par(oma = c(0,4,2,0))
dat<- radarchart(Buko_Aa, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Bukowno),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Buko_Ah, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Kato_Aa, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Katowice),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Kato_Ah, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Mias_Aa, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Miastezko~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Mias_Ah, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Piek_Aa, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Piekary~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Piek_Ah, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("red",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Zapa_Aa, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("black",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Zakopane),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Zapa_Ah, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("black",5)), pty=32, plty=c(3,2,1,2,3))

dat<- radarchart(Kowa_Aa, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("black",5)), pty=32, plty=c(3,2,1,2,3))
title(expression(Kowary),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(Kowa_Ah, axistype = 1, axislabcol="grey30", cglcol="grey40",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(rep("black",5)), pty=32, plty=c(3,2,1,2,3))

mtext(expression(italic(A.~arenosa)),side=2,line=1,outer=TRUE,cex=1.5,las=0,adj=0.8)
mtext(expression(italic(A.~halleri)),side=2,line=1,outer=TRUE,cex=1.5,las=0,adj=0.2)

dev.off()

pdf("X_Introgression_minall_species_separated_wKato2.pdf",width=20,height=8,paper="special",pointsize=14)
layout(matrix(1:12, ncol = 6))
par(mar = c(0, 0.5,1,0.5))
par(oma = c(0,4,2,0))

dat<- radarchart(rbind(Kowa_Aa[1:2,],Kowa_Ah[5,],Kowa_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Kowary),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Kowa_Ah[1:2,],Kowa_Aa[5,],Kowa_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Zapa_Aa[1:2,],Zapa_Ah[5,],Zapa_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Zakopane),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Zapa_Ah[1:2,],Zapa_Aa[5,],Zapa_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Kato_Aa[1:2,],Kato_Ah[5,],Kato_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Katowice),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Kato_Ah[1:2,],Kato_Aa[5,],Kato_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Buko_Aa[1:2,],Buko_Ah[5,],Buko_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Bukowno),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Buko_Ah[1:2,],Buko_Aa[5,],Buko_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Piek_Aa[1:2,],Piek_Ah[5,],Piek_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Piekary~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Piek_Ah[1:2,],Piek_Aa[5,],Piek_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dat<- radarchart(rbind(Mias_Aa[1:2,],Mias_Ah[5,],Mias_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))
title(expression(Miastezko~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Mias_Ah[1:2,],Mias_Aa[5,],Mias_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13],expression("H"^"+")))

dev.off()


#################################
#Leaf concentrations#
#################################

Leaf_col <- c( which( colnames(icp)=="Leaf_Pb" ), which( colnames(icp)=="Leaf_Cd" ), which( colnames(icp)=="Leaf_Zn" ),which( colnames(icp)=="Leaf_Mn" ),
             which( colnames(icp)=="Leaf_K" ),which( colnames(icp)=="Leaf_S" ),which( colnames(icp)=="Leaf_P" ),which( colnames(icp)=="Leaf_B" ),which( colnames(icp)=="Leaf_Cr" ),
            which( colnames(icp)=="Leaf_Al" ),which( colnames(icp)=="Leaf_Ca" ),which( colnames(icp)=="Leaf_Fe" ),which( colnames(icp)=="Leaf_Cu" ))

names(icp[Leaf_col])
Max<-log10(max(as.numeric(apply(icp[,Leaf_col],2,function(x) {x/min(x,na.rm=T)})),na.rm=T))



###################################################################################
#Buko arenosa
infoALeaf <- icp[icp$Population=="Buko"&icp$Species=="arenosa",1:2]
ArawLeaf <- icp[icp$Population=="Buko"&icp$Species=="arenosa",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_BukoLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_BukoLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Buko",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_BukoLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Buko_Aa <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Buko_Aa <- data.frame(t(Buko_Aa))
names(Buko_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")

#Buko halleri
infoALeaf <- icp[icp$Population=="Buko"&icp$Species=="halleri",1:2]
ArawLeaf <- icp[icp$Population=="Buko"&icp$Species=="halleri",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_BukoLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_BukoLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Buko",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_BukoLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Buko_Ah <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Buko_Ah <- data.frame(t(Buko_Ah))
names(Buko_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")


#Kato arenosa
infoALeaf <- icp[icp$Population=="Kato"&icp$Species=="arenosa",1:2]
ArawLeaf <- icp[icp$Population=="Kato"&icp$Species=="arenosa",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_KatoLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7), H=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_KatoLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Kato",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KatoLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Kato_Aa <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kato_Aa <- data.frame(t(Kato_Aa))
names(Kato_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")

#Kato halleri
infoALeaf <- icp[icp$Population=="Kato"&icp$Species=="halleri",1:2]
ArawLeaf <- icp[icp$Population=="Kato"&icp$Species=="halleri",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_KatoLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_KatoLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Kato",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KatoLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Kato_Ah <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kato_Ah <- data.frame(t(Kato_Ah))
names(Kato_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")



#Mias arenosa
infoALeaf <- icp[icp$Population=="Mias"&icp$Species=="arenosa",1:2]
ArawLeaf <- icp[icp$Population=="Mias"&icp$Species=="arenosa",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_MiasLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_MiasLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Mias",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_MiasLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Mias_Aa <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Mias_Aa <- data.frame(t(Mias_Aa))
names(Mias_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")

#Mias halleri
infoALeaf <- icp[icp$Population=="Mias"&icp$Species=="halleri",1:2]
ArawLeaf <- icp[icp$Population=="Mias"&icp$Species=="halleri",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_MiasLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_MiasLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Mias",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_MiasLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Mias_Ah <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Mias_Ah <- data.frame(t(Mias_Ah))
names(Mias_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")


#Piek arenosa
infoALeaf <- icp[icp$Population=="Piek"&icp$Species=="arenosa",1:2]
ArawLeaf <- icp[icp$Population=="Piek"&icp$Species=="arenosa",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_PiekLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_PiekLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Piek",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_PiekLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Piek_Aa <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Piek_Aa <- data.frame(t(Piek_Aa))
names(Piek_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")

#Piek halleri
infoALeaf <- icp[icp$Population=="Piek"&icp$Species=="halleri",1:2]
ArawLeaf <- icp[icp$Population=="Piek"&icp$Species=="halleri",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_PiekLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_PiekLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Piek",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_PiekLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Piek_Ah <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Piek_Ah <- data.frame(t(Piek_Ah))
names(Piek_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")



#Kowa arenosa
infoALeaf <- icp[icp$Population=="Kowa"&icp$Species=="arenosa",1:2]
ArawLeaf <- icp[icp$Population=="Kowa"&icp$Species=="arenosa",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_KowaLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_KowaLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Kowa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KowaLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Kowa_Aa <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kowa_Aa <- data.frame(t(Kowa_Aa))
names(Kowa_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")

#Kowa halleri
infoALeaf <- icp[icp$Population=="Kowa"&icp$Species=="halleri",1:2]
ArawLeaf <- icp[icp$Population=="Kowa"&icp$Species=="halleri",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_KowaLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_KowaLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Kowa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_KowaLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Kowa_Ah <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Kowa_Ah <- data.frame(t(Kowa_Ah))
names(Kowa_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")


#Zapa arenosa
infoALeaf <- icp[icp$Population=="Zapa"&icp$Species=="arenosa",1:2]
ArawLeaf <- icp[icp$Population=="Zapa"&icp$Species=="arenosa",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_ZapaLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7)
for ( i in 8:ncol(ArawLeaf)){
  
  stats_ZapaLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Zapa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_ZapaLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Zapa_Aa <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Zapa_Aa <- data.frame(t(Zapa_Aa))
names(Zapa_Aa) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")

#Zapa halleri
infoALeaf <- icp[icp$Population=="Zapa"&icp$Species=="halleri",1:2]
ArawLeaf <- icp[icp$Population=="Zapa"&icp$Species=="halleri",Leaf_col]    

#Leaf quantiles
#this part is only for saving quantiles statistics
stats_ZapaLeaf       <- data.frame(Pb=rep("NA",7), Cd=rep("NA",7),Zn=rep("NA",7),Mn=rep("NA",7),Ni=rep("NA",7),K=rep("NA",7),S=rep("NA",7),
                                P=rep("NA",7),B=rep("NA",7),Cr=rep("NA",7),Al=rep("NA",7),Mg=rep("NA",7),Ca=rep("NA",7),Fe=rep("NA",7),
                                Cu=rep("NA",7))
for ( i in 8:ncol(ArawLeaf)){
  
  stats_ZapaLeaf  [,i-7] <- quantile(ArawLeaf[infoALeaf$Population=="Zapa",i], probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=T)
}
row.names(stats_ZapaLeaf  )<- c("0","0.1","0.25","0.5","0.75","0.9","1")

length(ArawLeaf)
head(ArawLeaf)

As <- data.frame(Pb = ((ArawLeaf[,1]/min(icp[,Leaf_col][,1], na.rm=T))),   Cd = ((ArawLeaf[,2]/min(icp[,Leaf_col][,2], na.rm=T))),     Zn = ((ArawLeaf[,3]/min(icp[,Leaf_col][,3], na.rm=T))),
                 Mn = ((ArawLeaf[,4]/min(icp[,Leaf_col][,4], na.rm=T))),   K = ((ArawLeaf[,5]/min(icp[,Leaf_col][,5], na.rm=T))),
                 S = ((ArawLeaf[,6]/min(icp[,Leaf_col][,6], na.rm=T))),    P = ((ArawLeaf[,7]/min(icp[,Leaf_col][,7], na.rm=T))),      B = ((ArawLeaf[,8]/min(icp[,Leaf_col][,8], na.rm=T))),
                 Cr = ((ArawLeaf[,9]/min(icp[,Leaf_col][,9], na.rm=T))), Al = ((ArawLeaf[,10]/min(icp[,Leaf_col][,10], na.rm=T))),   Ca = ((ArawLeaf[,11]/min(icp[,Leaf_col][,11], na.rm=T))),
			Fe = ((ArawLeaf[,12]/min(icp[,Leaf_col][,12], na.rm=T))),   Cu = ((ArawLeaf[,13]/min(icp[,Leaf_col][,13], na.rm=T))))
A <- log10(As)#scaled by division with minimum and log10 transformed, also added H from pH, scaled in same way (scaled across N and NM Populationulation and both species)
A <- cbind(infoALeaf,A)

Zapa_Ah <- data.frame(max = rep(Max,13), min= c(rep(0,13)), M_halleri0= apply(A[,3:ncol(A)], 2, quantile,probs=c(0), na.rm=T), M_halleri10= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.1), na.rm=T), M_halleri50= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.5), na.rm=T), M_halleri90= apply(A[,3:ncol(A)], 2, quantile,probs=c(0.9), na.rm=T), M_halleri100= apply(A[,3:ncol(A)], 2, quantile,probs=c(1), na.rm=T))
Zapa_Ah <- data.frame(t(Zapa_Ah))
names(Zapa_Ah) <- c("Pb","Cd","Zn","Mn","K","S","P","B","Cr","Al","Ca","Fe","Cu")


pdf("Leaf_Introgression_minall_species_separated_wKato.pdf",width=20,height=8,paper="special")

layout(matrix(1:12, ncol = 6))
par(mar = c(0, 0.5,1,0.5))
par(oma = c(0,4,2,0))

dat<- radarchart(rbind(Kowa_Aa[1:2,],Kowa_Ah[5,],Kowa_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))
title(expression(Kowary),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Kowa_Ah[1:2,],Kowa_Aa[5,],Kowa_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))

dat<- radarchart(rbind(Zapa_Aa[1:2,],Zapa_Ah[5,],Zapa_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))
title(expression(Zakopane),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Zapa_Ah[1:2,],Zapa_Aa[5,],Zapa_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))

dat<- radarchart(rbind(Kato_Aa[1:2,],Kato_Ah[5,],Kato_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))
title(expression(Katowice),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Kato_Ah[1:2,],Kato_Aa[5,],Kato_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))

dat<- radarchart(rbind(Buko_Aa[1:2,],Buko_Ah[5,],Buko_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))
title(expression(Bukowno),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Buko_Ah[1:2,],Buko_Aa[5,],Buko_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))

dat<- radarchart(rbind(Piek_Aa[1:2,],Piek_Ah[5,],Piek_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))
title(expression(Piekary~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Piek_Ah[1:2,],Piek_Aa[5,],Piek_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))

dat<- radarchart(rbind(Mias_Aa[1:2,],Mias_Ah[5,],Mias_Aa[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkblue",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkgreen", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))
title(expression(Miastezko~Slaskie),cex.main=2.5,font.main=1,line=0)
dat<- radarchart(rbind(Mias_Ah[1:2,],Mias_Aa[5,],Mias_Ah[3:7,]), axistype = 1, axislabcol="grey40", cglcol="grey60",caxislabels=seq(0,5,1), cglty=1, cglwd=1, plwd = c(1,2,2,3,2,2), vlcex = 2, calcex = 2, pcol = c(NA,rep("darkgreen",5)), pty=32, plty=c(1,3,2,1,2,3),pfcol=c(scales::alpha("darkblue", 0.3),NA,NA,NA,NA,NA),vlabels=c(colnames(Kowa_Aa)[1:13]))

dev.off()

pdf("Spidercharts_legend.pdf",width=15,height=5,paper="special",pointsize=20)
plot(1,1)
legend("bottomright",lty=c(1,2,3),c("Median","10/90%ile","Minimum/Maximum"),horiz=T,lwd=3)
dev.off()

####
#Zn#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa"|icp_df$Species=="hybrid",]
icp_df$Species<-droplevels(icp_df$Species)
icp_df$Population<-droplevels(icp_df$Population)
boxplot(log10(icp_df$Leaf_Zn)~icp_df$Species*icp_df$Population,las=2,ylim=c(1,5),border=c("green","darkblue","darkcyan"),xlab="",ylab=expression(Log[10](Leaf~Zn)*"[ppm]"),cex.lab=1.3,cex.axis=1.2,na.action=na.exclude)

pdf("Leaf_ICP_introgression_boxplots.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa","hybrid"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df$Leaf_Zn)~icp_df$Species*icp_df$Population,las=2,ylim=c(1,5),border=c(rep(c("darkgreen","darkblue","darkcyan"),6)),names=c(rep("",18)),xlab="",ylab=expression(Log[10]*"(Leaf Zn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:2,4:5,7:11,13:17)],las=2,ylim=c(1,5),border=c("darkgreen","darkblue","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10]*"(Leaf Zn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
abline(h=log10(3000),col="red")
stripchart(log10(icp_df$Leaf_Zn)~icp_df$Species*icp_df$Population,at=c(1:2,NA,NA,4,NA,5:9,NA,10:12,13:14,NA),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue","darkcyan"),6)))
stripchart(log10(icp_df$Leaf_Zn)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df$Leaf_Zn)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,0.75, 2.4,0.75, xpd = TRUE,lwd=2)
segments(2.6,0.75, 4.4,0.75, xpd = TRUE,lwd=2)
segments(4.6,0.75, 7.4,0.75, xpd = TRUE,lwd=2)
segments(7.6,0.75, 9.4,0.75, xpd = TRUE,lwd=2)
segments(9.6,0.75, 12.4,0.75, xpd = TRUE,lwd=2)
segments(12.6,0.75, 14.4,0.75, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,6,8.5,11,13.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue","darkcyan"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa),italic(Hybrid))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df2$Leaf_Zn)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(1,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Zn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(1,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Zn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
abline(h=log10(3000),col="red")
stripchart(log10(icp_df2$Leaf_Zn)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,0.75, 2.4,0.75, xpd = TRUE,lwd=2)
segments(2.6,0.75, 4.4,0.75, xpd = TRUE,lwd=2)
segments(4.6,0.75, 6.4,0.75, xpd = TRUE,lwd=2)
segments(6.6,0.75, 8.4,0.75, xpd = TRUE,lwd=2)
segments(8.6,0.75, 10.4,0.75, xpd = TRUE,lwd=2)
segments(10.6,0.75, 12.4,0.75, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

require(Hmisc)

pdf("Leaf_ICP_introgression_boxplots_wohybrids_Zn_neworder.pdf",width=11,height=10.2,paper="special",pointsize=25)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Kowa","Zapa","Kato","Buko","Piek","Mias"))
par(lwd=2.5)
par(mar=c(5,6.5,1,1))
par(mgp=c(4.5,0.75,0))
par(oma=c(1,1,1,1))

x<-boxplot(log10(icp_df2$Leaf_Zn)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(1,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Leaf~Zn~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5,plot=F)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(1,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Leaf~Zn~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n",xlim=c(0.9,12.1))
abline(h=log10(3000),col="grey40",lwd=2.5)
stripchart(log10(icp_df2$Leaf_Zn)~icp_df2$Species*icp_df2$Population,at=c(1:4,NA,6:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))
mgp.axis(side=2,at=c(1:5),labels=c("10","100","1,000","10,000","100,000"),las=1,cex.axis=1.3)

segments(0.6,0.75, 2.4,0.75, xpd = TRUE,lwd=2.5)
segments(2.6,0.75, 4.4,0.75, xpd = TRUE,lwd=2.5)
segments(4.6,0.75, 6.4,0.75, xpd = TRUE,lwd=2.5)
segments(6.6,0.75, 8.4,0.75, xpd = TRUE,lwd=2.5)
segments(8.6,0.75, 10.4,0.75, xpd = TRUE,lwd=2.5)
segments(10.6,0.75, 12.4,0.75, xpd = TRUE,lwd=2.5)

mtext(c("Kowa","Zapa","Kato","Buko","Piek","Mias"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c("black","black",rep("red",4)),font=2)
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=6,horiz=T,xpd=T,inset=c(-0.25,-0.3),text.width=4.5)
dev.off()


####
#Cd#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa"|icp_df$Species=="hybrid",]
icp_df$Species<-droplevels(icp_df$Species)
icp_df$Population<-droplevels(icp_df$Population)

pdf("Leaf_ICP_introgression_boxplots_Cd.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa","hybrid"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))
x<-boxplot(log10(icp_df$Leaf_Cd+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,4),border=c(rep(c("darkgreen","darkblue","darkcyan"),6)),names=c(rep("",18)),xlab="",ylab=expression(Log[10](Leaf~Cd+1)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
boxplot(x$stats[,c(1:2,4:5,7:11,13:17)],las=2,ylim=c(0,4),border=c("darkgreen","darkblue","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Cd+1)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
abline(h=log10(100+1),col="red")
stripchart(log10(icp_df$Leaf_Cd+1)~icp_df$Species*icp_df$Population,at=c(1:2,NA,3:4,NA,5:9,NA,10:12,13:14,NA),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE)
segments(4.6,-0.25, 7.4,-0.25, xpd = TRUE)
segments(7.6,-0.25, 9.4,-0.25, xpd = TRUE)
segments(9.6,-0.25, 12.4,-0.25, xpd = TRUE)
segments(12.6,-0.25, 14.4,-0.25, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,6,8.5,11,13.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue","darkcyan"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa),italic(Hybrid))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_Cd_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df2$Leaf_Cd+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Cd [ppm] +1)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Cd [ppm] +1)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
abline(h=log10(100+1),col="red")
stripchart(log10(icp_df2$Leaf_Cd+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Cd+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Cd+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,-0.2, 2.4,-0.2, xpd = TRUE,lwd=2)
segments(2.6,-0.2, 4.4,-0.2, xpd = TRUE,lwd=2)
segments(4.6,-0.2, 6.4,-0.2, xpd = TRUE,lwd=2)
segments(6.6,-0.2, 8.4,-0.2, xpd = TRUE,lwd=2)
segments(8.6,-0.2, 10.4,-0.2, xpd = TRUE,lwd=2)
segments(10.6,-0.2, 12.4,-0.2, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()


pdf("Leaf_ICP_introgression_boxplots_Cd_neworder.pdf",width=11,height=10.2,paper="special",pointsize=25)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Kowa","Zapa","Kato","Buko","Piek","Mias"))
par(lwd=2.5)
par(mar=c(5,6.5,1,1))
par(mgp=c(4,0.75,0))
par(oma=c(1,1,1,1))
x<-boxplot(log10(icp_df2$Leaf_Cd+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,4),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab=expression(Leaf~Cd~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,plot=F)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Leaf~Cd~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n",xlim=c(0.9,12.1))
abline(h=log10(100+1),col="grey40")
stripchart(log10(icp_df2$Leaf_Cd+1)~icp_df2$Species*icp_df2$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Cd+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Cd+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))
mgp.axis(side=2,at=c(1:4),labels=c("10","100","1,000","10,000"),las=1,cex.axis=1.3)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE,lwd=2.5)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE,lwd=2.5)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE,lwd=2.5)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE,lwd=2.5)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE,lwd=2.5)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE,lwd=2.5)

mtext(c("Kowa","Zapa","Kato","Buko","Piek","Mias"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c("black","black",rep("red",4)),font=2)
par(lwd=2.5)
#legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

####
#Pb#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa"|icp_df$Species=="hybrid",]
icp_df$Species<-droplevels(icp_df$Species)
icp_df$Population<-droplevels(icp_df$Population)

pdf("Leaf_ICP_introgression_boxplots_Pb.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa","hybrid"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))
x<-boxplot(log10(icp_df$Leaf_Pb+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,4),border=c(rep(c("darkgreen","darkblue","darkcyan"),6)),names=c(rep("",18)),xlab="",ylab=expression(Log[10](Leaf~Pb+1)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
boxplot(x$stats[,c(1:2,4:5,7:11,13:17)],las=2,ylim=c(0,4),border=c("darkgreen","darkblue","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10](Leaf~Pb+1)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
abline(h=log10(1000+1),col="red")
stripchart(log10(icp_df$Leaf_Pb+1)~icp_df$Species*icp_df$Population,at=c(1:2,NA,3:4,NA,5:9,NA,10:12,13:14,NA),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE)
segments(4.6,-0.25, 7.4,-0.25, xpd = TRUE)
segments(7.6,-0.25, 9.4,-0.25, xpd = TRUE)
segments(9.6,-0.25, 12.4,-0.25, xpd = TRUE)
segments(12.6,-0.25, 14.4,-0.25, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,6,8.5,11,13.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue","darkcyan"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa),italic(Hybrid))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_Pb_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df2$Leaf_Pb+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Pb [ppm] +1)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Pb [ppm] +1)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
abline(h=log10(100+1),col="red")
stripchart(log10(icp_df2$Leaf_Pb+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Pb+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Pb+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE,lwd=2)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE,lwd=2)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE,lwd=2)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE,lwd=2)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE,lwd=2)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()


pdf("Leaf_ICP_introgression_boxplots_Pb_neworder.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Kowa","Zapa","Kato","Buko","Piek","Mias"))
par(lwd=1)
par(mar=c(6,5,1,1))
x<-boxplot(log10(icp_df2$Leaf_Pb+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,4),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab=expression(Log[10](Leaf~Pb+1)*"[ppm]"),xaxt="n",cex.lab=1.5,cex.axis=1.3)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Pb [ppm] +1)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
abline(h=log10(1000+1),col="red")
stripchart(log10(icp_df2$Leaf_Pb+1)~icp_df2$Species*icp_df2$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Pb+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Pb+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE,lwd=2)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE,lwd=2)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE,lwd=2)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE,lwd=2)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE,lwd=2)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE,lwd=2)

mtext(c("Kowa","Zapa","Kato","Buko","Piek","Mias"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c("black","black",rep("red",4)))
par(lwd=2)
#legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()


##################################
#Ratios#
##################################

####
#Zn#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa",]
icp_df$Species<-droplevels(icp_df$Species)
icp_df$Population<-droplevels(icp_df$Population)

pdf("Leaf_ICP_introgression_boxplots_Zn_leafbyextractableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Kowa","Zapa","Kato","Buko","Piek","Mias"))
par(lwd=2.5)
par(mar=c(5,6,1,1))
par(mgp=c(4,0.75,0))
par(oma=c(1,1,1,1))

x<-boxplot(log10(icp_df2$Leaf_Zn/icp_df$X_Zn)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(-2,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression("Leaf Zn conc. (ppm)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5,plot=F)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(-2,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n",xlim=c(0.9,12.1))
stripchart(log10(icp_df2$Leaf_Zn/icp_df$X_Zn)~icp_df2$Species*icp_df2$Population,at=c(1:4,NA,6:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Zn/icp_df$X_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Zn/icp_df$X_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))
axis(side=2,at=c(-2:5),labels=c("0.01","0.1","1","10","100","1,000","10,000","100,000"),las=1,cex.axis=1.3)

mtext("(Zn conc. in leaves/ Zn conc. in soil)",side=2,line=3,cex=1)
mtext("Zn accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-2.3, 2.4,-2.3, xpd = TRUE,lwd=2.5)
segments(2.6,-2.3, 4.4,-2.3, xpd = TRUE,lwd=2.5)
segments(4.6,-2.3, 6.4,-2.3, xpd = TRUE,lwd=2.5)
segments(6.6,-2.3, 8.4,-2.3, xpd = TRUE,lwd=2.5)
segments(8.6,-2.3, 10.4,-2.3, xpd = TRUE,lwd=2.5)
segments(10.6,-2.3, 12.4,-2.3, xpd = TRUE,lwd=2.5)

mtext(c("Kowa","Zapa","Kato","Buko","Piek","Mias"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c("black","black",rep("red",4)),font=2)
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4,horiz=T,xpd=T,inset=c(-0.25,-0.3),text.width=4.5)
dev.off()


icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa",]
icp_df$Species<-droplevels(icp_df$Species)
icp_df$Population<-droplevels(icp_df$Population)

pdf("Leaf_ICP_introgression_boxplots_Zn_leafbyexchangeableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Kowa","Zapa","Kato","Buko","Piek","Mias"))
par(lwd=2.5)
par(mar=c(5,6,1,1))
par(mgp=c(4,0.75,0))
par(oma=c(1,1,1,1))

x<-boxplot(log10(icp_df2$Leaf_Zn/icp_df$CH_Zn)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression("Leaf Zn conc. (ppm)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5,plot=F)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n",xlim=c(0.9,12.1))
stripchart(log10(icp_df2$Leaf_Zn/icp_df$CH_Zn)~icp_df2$Species*icp_df2$Population,at=c(1:4,NA,6:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Zn/icp_df$CH_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Zn/icp_df$CH_Zn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))
axis(side=2,at=c(0:5),labels=c("0.1","1","10","100","1,000","10,000"),las=1,cex.axis=1.3)

mtext("(Zn conc. in leaves/ Zn conc. in soil)",side=2,line=3,cex=1)
mtext("Zn accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.3, 2.4,-0.3, xpd = TRUE,lwd=2.5)
segments(2.6,-0.3, 4.4,-0.3, xpd = TRUE,lwd=2.5)
segments(4.6,-0.3, 6.4,-0.3, xpd = TRUE,lwd=2.5)
segments(6.6,-0.3, 8.4,-0.3, xpd = TRUE,lwd=2.5)
segments(8.6,-0.3, 10.4,-0.3, xpd = TRUE,lwd=2.5)
segments(10.6,-0.3, 12.4,-0.3, xpd = TRUE,lwd=2.5)

mtext(c("Kowa","Zapa","Kato","Buko","Piek","Mias"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c("black","black",rep("red",4)),font=2)
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4,horiz=T,xpd=T,inset=c(-0.25,-0.3),text.width=4.5)
dev.off()



icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))
x<-boxplot(log10((icp_df$Leaf_Zn/icp_df$X_Zn)+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,4),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
boxplot(x$stats,las=2,ylim=c(0,4),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
stripchart(log10((icp_df$Leaf_Zn/icp_df$X_Zn)+1)~icp_df$Species*icp_df$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

mtext("(Zn conc. in leaves/ Zn conc. in soil)",side=2,line=3,cex=1)
mtext("Zn accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)),border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_Zn_leafbyextractableratio_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))

x<-boxplot(log10((icp_df2$Leaf_Zn/icp_df2$X_Zn)+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,4),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10((icp_df2$Leaf_Zn/icp_df2$X_Zn)+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10((icp_df2$Leaf_Zn/icp_df2$X_Zn)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10((icp_df2$Leaf_Zn/icp_df2$X_Zn)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

mtext("(Zn conc. in leaves/ Zn conc. in soil)",side=2,line=3,cex=1)
mtext("Zn accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE,lwd=2)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE,lwd=2)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE,lwd=2)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE,lwd=2)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE,lwd=2)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

####
#Cd#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa",]
icp_df$Species<-droplevels(icp_df$Species)
icp_df$Population<-droplevels(icp_df$Population)

pdf("Leaf_ICP_introgression_boxplots_Cd_leafbyextractableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Kowa","Zapa","Kato","Buko","Piek","Mias"))
par(lwd=2.5)
par(mar=c(5,6,1,1))
par(mgp=c(4,0.75,0))
par(oma=c(1,1,1,1))

x<-boxplot(log10(icp_df2$Leaf_Cd/icp_df$X_Cd)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(-2,3),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression("Leaf Cd conc. (ppm)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5,plot=F)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(-2,3),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n",xlim=c(0.9,12.1))
stripchart(log10(icp_df2$Leaf_Cd/icp_df$X_Cd)~icp_df2$Species*icp_df2$Population,at=c(1:4,NA,6:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Cd/icp_df$X_Cd)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Cd/icp_df$X_Cd)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))
axis(side=2,at=c(-2:4),labels=c("0.01","0.1","1","10","100","1,000","10,000"),las=1,cex.axis=1.3)

mtext("(Cd conc. in leaves/ Cd conc. in soil)",side=2,line=3,cex=1)
mtext("Cd accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-2.3, 2.4,-2.3, xpd = TRUE,lwd=2.5)
segments(2.6,-2.3, 4.4,-2.3, xpd = TRUE,lwd=2.5)
segments(4.6,-2.3, 6.4,-2.3, xpd = TRUE,lwd=2.5)
segments(6.6,-2.3, 8.4,-2.3, xpd = TRUE,lwd=2.5)
segments(8.6,-2.3, 10.4,-2.3, xpd = TRUE,lwd=2.5)
segments(10.6,-2.3, 12.4,-2.3, xpd = TRUE,lwd=2.5)

mtext(c("Kowa","Zapa","Kato","Buko","Piek","Mias"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c("black","black",rep("red",4)),font=2)
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4,horiz=T,xpd=T,inset=c(-0.25,-0.3),text.width=4.5)
dev.off()


pdf("Leaf_ICP_introgression_boxplots_Cd_leafbyexchangeableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Kowa","Zapa","Kato","Buko","Piek","Mias"))
par(lwd=2.5)
par(mar=c(5,6,1,1))
par(mgp=c(4,0.75,0))
par(oma=c(1,1,1,1))

x<-boxplot(log10(icp_df2$Leaf_Cd/icp_df$CH_Cd)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression("Leaf Cd conc. (ppm)"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5,plot=F)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n",xlim=c(0.9,12.1))
stripchart(log10(icp_df2$Leaf_Cd/icp_df$CH_Cd)~icp_df2$Species*icp_df2$Population,at=c(1:4,NA,6:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Cd/icp_df$CH_Cd)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Cd/icp_df$CH_Cd)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(5),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))
axis(side=2,at=c(0:5),labels=c("0.1","1","10","100","1,000","10,000"),las=1,cex.axis=1.3)

mtext("(Cd conc. in leaves/ Cd conc. in soil)",side=2,line=3,cex=1)
mtext("Cd accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.3, 2.4,-0.3, xpd = TRUE,lwd=2.5)
segments(2.6,-0.3, 4.4,-0.3, xpd = TRUE,lwd=2.5)
segments(4.6,-0.3, 6.4,-0.3, xpd = TRUE,lwd=2.5)
segments(6.6,-0.3, 8.4,-0.3, xpd = TRUE,lwd=2.5)
segments(8.6,-0.3, 10.4,-0.3, xpd = TRUE,lwd=2.5)
segments(10.6,-0.3, 12.4,-0.3, xpd = TRUE,lwd=2.5)

mtext(c("Kowa","Zapa","Kato","Buko","Piek","Mias"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c("black","black",rep("red",4)),font=2)
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4,horiz=T,xpd=T,inset=c(-0.25,-0.3),text.width=4.5)
dev.off()




pdf("Leaf_ICP_introgression_boxplots_Cd_leafbyextractableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))
x<-boxplot(log10((icp_df$Leaf_Cd/icp_df$X_Cd)+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,3),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
boxplot(x$stats,las=2,ylim=c(0,3),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
stripchart(log10((icp_df$Leaf_Cd/icp_df$X_Cd)+1)~icp_df$Species*icp_df$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

mtext("(Cd conc. in leaves/ Cd conc. in soil)",side=2,line=3,cex=1)
mtext("Cd accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.2, 2.4,-0.2, xpd = TRUE)
segments(2.6,-0.2, 4.4,-0.2, xpd = TRUE)
segments(4.6,-0.2, 6.4,-0.2, xpd = TRUE)
segments(6.6,-0.2, 8.4,-0.2, xpd = TRUE)
segments(8.6,-0.2, 10.4,-0.2, xpd = TRUE)
segments(10.6,-0.2, 12.4,-0.2, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)),border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_Cd_leafbyextractableratio_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))

x<-boxplot(log10((icp_df2$Leaf_Cd/icp_df2$X_Cd)+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,3),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,3),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10((icp_df2$Leaf_Cd/icp_df2$X_Cd)+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10((icp_df2$Leaf_Cd/icp_df2$X_Cd)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10((icp_df2$Leaf_Cd/icp_df2$X_Cd)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

mtext("(Cd conc. in leaves/ Cd conc. in soil)",side=2,line=3,cex=1)
mtext("Cd accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.2, 2.4,-0.2, xpd = TRUE,lwd=2)
segments(2.6,-0.2, 4.4,-0.2, xpd = TRUE,lwd=2)
segments(4.6,-0.2, 6.4,-0.2, xpd = TRUE,lwd=2)
segments(6.6,-0.2, 8.4,-0.2, xpd = TRUE,lwd=2)
segments(8.6,-0.2, 10.4,-0.2, xpd = TRUE,lwd=2)
segments(10.6,-0.2, 12.4,-0.2, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

####
#Pb#
####
pdf("Leaf_ICP_introgression_boxplots_Pb_leafbyextractableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))
x<-boxplot(log10((icp_df$Leaf_Pb/icp_df$X_Pb)+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,2.5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
boxplot(x$stats,las=2,ylim=c(0,2.5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
stripchart(log10((icp_df$Leaf_Pb/icp_df$X_Pb)+1)~icp_df$Species*icp_df$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

mtext("(Pb conc. in leaves/ Pb conc. in soil)",side=2,line=3,cex=1)
mtext("Pb accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.15, 2.4,-0.15, xpd = TRUE)
segments(2.6,-0.15, 4.4,-0.15, xpd = TRUE)
segments(4.6,-0.15, 6.4,-0.15, xpd = TRUE)
segments(6.6,-0.15, 8.4,-0.15, xpd = TRUE)
segments(8.6,-0.15, 10.4,-0.15, xpd = TRUE)
segments(10.6,-0.15, 12.4,-0.15, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)),border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()


pdf("Leaf_ICP_introgression_boxplots_Pb_leafbyextractableratio_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))

x<-boxplot(log10((icp_df2$Leaf_Pb/icp_df2$X_Pb)+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,2.5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,2.5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10((icp_df2$Leaf_Pb/icp_df2$X_Pb)+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10((icp_df2$Leaf_Pb/icp_df2$X_Pb)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10((icp_df2$Leaf_Pb/icp_df2$X_Pb)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

mtext("(Pb conc. in leaves/ Pb conc. in soil)",side=2,line=3,cex=1)
mtext("Pb accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.15, 2.4,-0.15, xpd = TRUE,lwd=2)
segments(2.6,-0.15, 4.4,-0.15, xpd = TRUE,lwd=2)
segments(4.6,-0.15, 6.4,-0.15, xpd = TRUE,lwd=2)
segments(6.6,-0.15, 8.4,-0.15, xpd = TRUE,lwd=2)
segments(8.6,-0.15, 10.4,-0.15, xpd = TRUE,lwd=2)
segments(10.6,-0.15, 12.4,-0.15, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

##################################
#Ratios exchangeable soil conc.#
##################################

####
#Zn#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa",]
icp_df$Species<-droplevels(icp_df$Species)
icp_df$Population<-droplevels(icp_df$Population)

pdf("Leaf_ICP_introgression_boxplots_Zn_leafbyexchangeableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))
x<-boxplot(log10((icp_df$Leaf_Zn/icp_df$CH_Zn)+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
boxplot(x$stats,las=2,ylim=c(0,5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
stripchart(log10((icp_df$Leaf_Zn/icp_df$CH_Zn)+1)~icp_df$Species*icp_df$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

mtext("(Zn conc. in leaves/ Zn conc. in soil)",side=2,line=3,cex=1)
mtext("Zn accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)),border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_Zn_leafbyexchangeableratio_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))

x<-boxplot(log10((icp_df2$Leaf_Zn/icp_df2$CH_Zn)+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10((icp_df2$Leaf_Zn/icp_df2$CH_Zn)+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10((icp_df2$Leaf_Zn/icp_df2$CH_Zn)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10((icp_df2$Leaf_Zn/icp_df2$CH_Zn)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

mtext("(Zn conc. in leaves/ Zn conc. in soil)",side=2,line=3,cex=1)
mtext("Zn accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE,lwd=2)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE,lwd=2)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE,lwd=2)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE,lwd=2)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE,lwd=2)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

####
#Cd#
####
pdf("Leaf_ICP_introgression_boxplots_Cd_leafbyexchangeableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))
x<-boxplot(log10((icp_df$Leaf_Cd/icp_df$CH_Cd)+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
boxplot(x$stats,las=2,ylim=c(0,5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
stripchart(log10((icp_df$Leaf_Cd/icp_df$CH_Cd)+1)~icp_df$Species*icp_df$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

mtext("(Cd conc. in leaves/ Cd conc. in soil)",side=2,line=3,cex=1)
mtext("Cd accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)),border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_Cd_leafbyexchangeableratio_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))

x<-boxplot(log10((icp_df2$Leaf_Cd/icp_df2$CH_Cd)+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10((icp_df2$Leaf_Cd/icp_df2$CH_Cd)+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10((icp_df2$Leaf_Cd/icp_df2$CH_Cd)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10((icp_df2$Leaf_Cd/icp_df2$CH_Cd)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

mtext("(Cd conc. in leaves/ Cd conc. in soil)",side=2,line=3,cex=1)
mtext("Cd accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE,lwd=2)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE,lwd=2)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE,lwd=2)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE,lwd=2)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE,lwd=2)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()


####
#Pb#
####
pdf("Leaf_ICP_introgression_boxplots_Pb_leafbyexchangeableratio.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))
x<-boxplot(log10((icp_df$Leaf_Pb/icp_df$CH_Pb)+1)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
boxplot(x$stats,las=2,ylim=c(0,5),border=c(rep(c("darkgreen","darkblue"),6)),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)))
stripchart(log10((icp_df$Leaf_Pb/icp_df$CH_Pb)+1)~icp_df$Species*icp_df$Population,at=c(1:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col='black')

mtext("(Pb conc. in leaves/ Pb conc. in soil)",side=2,line=3,cex=1)
mtext("Pb accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill=c(rgb(0,100,0,max=255,alpha=130),rgb(0,0,139,max=255,alpha=130)),border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_Pb_leafbyexchangeableratio_wohybrids.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,6,1,1))

x<-boxplot(log10((icp_df2$Leaf_Pb/icp_df2$CH_Pb)+1)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab="",xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10((icp_df2$Leaf_Pb/icp_df2$CH_Pb)+1)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10((icp_df2$Leaf_Pb/icp_df2$CH_Pb)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10((icp_df2$Leaf_Pb/icp_df2$CH_Pb)+1)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

mtext("(Pb conc. in leaves/ Pb conc. in soil)",side=2,line=3,cex=1)
mtext("Pb accumulation efficiency",side=2,line=4,cex=1.5)

segments(0.6,-0.25, 2.4,-0.25, xpd = TRUE,lwd=2)
segments(2.6,-0.25, 4.4,-0.25, xpd = TRUE,lwd=2)
segments(4.6,-0.25, 6.4,-0.25, xpd = TRUE,lwd=2)
segments(6.6,-0.25, 8.4,-0.25, xpd = TRUE,lwd=2)
segments(8.6,-0.25, 10.4,-0.25, xpd = TRUE,lwd=2)
segments(10.6,-0.25, 12.4,-0.25, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()


#########################################
#Leaf against soil
#########################################

icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))

plot(log10(icp_df$Leaf_Zn)~log10(icp_df$X_Zn),pch=c(21:25,16)[icp_df$Population],col=c("darkgreen","darkblue")[icp_df$Species],cex=1.5,bg=c("darkgreen","darkblue")[icp_df$Species],ylab=expression(Log[10](Leaf~Zn)*"[ppm]"),xlab=expression(Log[10](Soil~Zn)*"[ppm]"))
legend(0,4.7,c("Buko","Kato","Mias","Piek","Zapa","Kowa"),pch=c(21:25,16))
legend(0,3.8,c(expression(italic(A.~halleri),italic(A.~arenosa))),fill=c("darkgreen","darkblue"))

icp_dfA<-icp_df[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID!="Kato_h09",]
icp_dfB<-icp_df[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID=="Kato_h09",]
icp_dfC<-icp_df[!(icp_df$Species=="halleri"&icp_df$Population=="Kato"),]

pdf("ICP_introgression_Zn_leaf_X_wKato.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_dfC$Leaf_Zn)~log10(icp_dfC$X_Zn),pch=c(15,17)[icp_dfC$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_dfC$Population],cex=1.5,ylab=expression(Log[10]*"(Leaf Zn [ppm])"),xlab=expression(Log[10]*"(Soil Zn [ppm])"),cex.lab=1.5,cex.axis=1.3)
points(log10(icp_dfA$Leaf_Zn)~log10(icp_dfA$X_Zn),pch=22,col="darkblue",bg="white",cex=1.5,lwd=2)
points(log10(icp_dfB$Leaf_Zn)~log10(icp_dfB$X_Zn),pch=15,col="darkblue",cex=1.5)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()

pdf("ICP_introgression_Cd_leaf_X.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_df$Leaf_Cd+1)~log10(icp_df$X_Cd+1),pch=c(15,17)[icp_df$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_df$Population],cex=1.5,bg=c("darkgreen","darkblue")[icp_df$Species],ylab=expression(Log[10]*"(Leaf Cd [ppm] +1)"),xlab=expression(Log[10]*"(Soil Cd [ppm] +1)"),cex.lab=1.5,cex.axis=1.3)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()

pdf("ICP_introgression_Cd_leaf_X_wKato.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_dfC$Leaf_Cd)~log10(icp_dfC$X_Cd),pch=c(15,17)[icp_dfC$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_dfC$Population],cex=1.5,ylab=expression(Log[10]*"(Leaf Cd [ppm])"),xlab=expression(Log[10]*"(Soil Cd [ppm])"),cex.lab=1.5,cex.axis=1.3)
points(log10(icp_dfA$Leaf_Cd)~log10(icp_dfA$X_Cd),pch=22,col="darkblue",bg="white",cex=1.5,lwd=2)
points(log10(icp_dfB$Leaf_Cd)~log10(icp_dfB$X_Cd),pch=15,col="darkblue",cex=1.5)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()

pdf("ICP_introgression_Pb_leaf_X.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_df$Leaf_Pb+1)~log10(icp_df$X_Pb+1),pch=c(15,17)[icp_df$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_df$Population],cex=1.5,bg=c("darkgreen","darkblue")[icp_df$Species],ylab=expression(Log[10]*"(Leaf Pb [ppm] +1)"),xlab=expression(Log[10]*"(Soil Pb [ppm] +1)"),cex.lab=1.5,cex.axis=1.3)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()

pdf("ICP_introgression_Pb_leaf_X_wKato.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_dfC$Leaf_Pb)~log10(icp_dfC$X_Pb),pch=c(15,17)[icp_dfC$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_dfC$Population],cex=1.5,ylab=expression(Log[10]*"(Leaf Pb [ppm])"),xlab=expression(Log[10]*"(Soil Pb [ppm])"),cex.lab=1.5,cex.axis=1.3)
points(log10(icp_dfA$Leaf_Pb)~log10(icp_dfA$X_Pb),pch=22,col="darkblue",bg="white",cex=1.5,lwd=2)
points(log10(icp_dfB$Leaf_Pb)~log10(icp_dfB$X_Pb),pch=15,col="darkblue",cex=1.5)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()



pdf("ICP_introgression_Zn_leaf_CH.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,1))
plot(log10(icp_df$Leaf_Zn+1)~log10(icp_df$CH_Zn+1),pch=c(15,17)[icp_df$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_df$Population],cex=1.5,bg=c("darkgreen","darkblue")[icp_df$Species],ylab=expression(Log[10]*"(Leaf Zn [ppm] +1)"),xlab=expression(Log[10]*"(Soil Zn [ppm] +1)"),cex.lab=1.5,cex.axis=1.3)
legend(2.1,3,c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4)
legend(2.1,1.95,c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("ICP_introgression_Zn_leaf_CH_wKato.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_dfC$Leaf_Zn)~log10(icp_dfC$CH_Zn),pch=c(15,17)[icp_dfC$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_dfC$Population],cex=1.5,ylab=expression(Log[10]*"(Leaf Zn [ppm])"),xlab=expression(Log[10]*"(Soil Zn [ppm])"),cex.lab=1.5,cex.axis=1.3)
points(log10(icp_dfA$Leaf_Zn)~log10(icp_dfA$CH_Zn),pch=22,col="darkblue",bg="white",cex=1.5,lwd=2)
points(log10(icp_dfB$Leaf_Zn)~log10(icp_dfB$CH_Zn),pch=15,col="darkblue",cex=1.5)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()


pdf("ICP_introgression_Cd_leaf_CH.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,1))
plot(log10(icp_df$Leaf_Cd+1)~log10(icp_df$CH_Cd+1),pch=c(15,17)[icp_df$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_df$Population],cex=1.5,bg=c("darkgreen","darkblue")[icp_df$Species],ylab=expression(Log[10]*"(Leaf Cd [ppm] +1)"),xlab=expression(Log[10]*"(Soil Cd [ppm] +1)"),cex.lab=1.5,cex.axis=1.3)
legend(1.7,1.6,c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4)
legend(1.7,0.5,c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("ICP_introgression_Cd_leaf_CH_wKato.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_dfC$Leaf_Cd)~log10(icp_dfC$CH_Cd),pch=c(15,17)[icp_dfC$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_dfC$Population],cex=1.5,ylab=expression(Log[10]*"(Leaf Cd [ppm])"),xlab=expression(Log[10]*"(Soil Cd [ppm])"),cex.lab=1.5,cex.axis=1.3)
points(log10(icp_dfA$Leaf_Cd)~log10(icp_dfA$CH_Cd),pch=22,col="darkblue",bg="white",cex=1.5,lwd=2)
points(log10(icp_dfB$Leaf_Cd)~log10(icp_dfB$CH_Cd),pch=15,col="darkblue",cex=1.5)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()


pdf("ICP_introgression_Pb_leaf_CH.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,1))
plot(log10(icp_df$Leaf_Pb+1)~log10(icp_df$CH_Pb+1),pch=c(15,17)[icp_df$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_df$Population],cex=1.5,bg=c("darkgreen","darkblue")[icp_df$Species],ylab=expression(Log[10]*"(Leaf Pb [ppm] +1)"),xlab=expression(Log[10]*"(Soil Pb [ppm] +1)"),cex.lab=1.5,cex.axis=1.3)
legend(1.2,1.75,c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4)
legend(1.2,0.45,c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("ICP_introgression_Pb_leaf_CH_wKato.pdf",width=8,height=8,paper="special",pointsize=14)
par(lwd=1)
par(mar=c(6,6,1,10))
plot(log10(icp_dfC$Leaf_Pb)~log10(icp_dfC$CH_Pb),pch=c(15,17)[icp_dfC$Species],col=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue")[icp_dfC$Population],cex=1.5,ylab=expression(Log[10]*"(Leaf Pb [ppm])"),xlab=expression(Log[10]*"(Soil Pb [ppm])"),cex.lab=1.5,cex.axis=1.3)
points(log10(icp_dfA$Leaf_Pb)~log10(icp_dfA$CH_Pb),pch=22,col="darkblue",bg="white",cex=1.5,lwd=2)
points(log10(icp_dfB$Leaf_Pb)~log10(icp_dfB$CH_Pb),pch=15,col="darkblue",cex=1.5)
legend("bottomright",c("Buko","Kato","Mias","Piek","Zapa","Kowa"),fill=c("darkgreen","darkblue","darkred","orange","lightgreen","lightblue"),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.35,0.2),xpd=T)
legend("bottomright",c(expression(italic(A.~halleri),italic(A.~arenosa))),pch=c(15,17),cex=1.5,bty="n",pt.lwd=4,inset=c(-0.5,0),xpd=T)
dev.off()



#########################################################
#Subsets identification#
#########################################################

mean(icp$Leaf_Zn[icp$Species=="arenosa"&icp$Population=="Kato"])
#326.2848
#*3=978.8545

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Buko")&(icp$Leaf_Zn>3000))]
#Buko_a28

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Buko")&(icp$Leaf_Zn<3000)&(icp$Leaf_Zn>1000))]
#"Buko_a22"  "Buko_a30"  "Buko_a34"  "Buko_a34b"

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Kato")&(icp$Leaf_Zn>3000))]
#0

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Piek")&(icp$Leaf_Zn>3000))]
# "Piek_a1"  "Piek_a10" "Piek_a12" "Piek_a13" "Piek_a14" "Piek_a15" "Piek_a2" "Piek_a4"  "Piek_a5"  "Piek_a7"  "Piek_a8" "Piek_a9"

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Mias")&(icp$Leaf_Zn>3000))]
"Mias_a43"   "Mias_a44"  "Mias_a45"   "Mias_a46"   "Mias_a47"   "Mias_a48"   "Mias_a50"  "Mias_a55"   "Mias001a18" "Mias003a09" "Mias003a10" "Mias003a11""Mias003a13" "Mias003a15" "Mias003a16"

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Piek")&(icp$Leaf_Zn<3000))]
"Piek_a16"   "Piek_a17"   "Piek_a18"   "Piek_a19"   "Piek_a20-1"	"Piek_a22"

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Zapa")&(icp$Leaf_Zn>3000))]
#0

icp$Sample_ID[((icp$Species=="arenosa")&(icp$Population=="Mias")&(icp$Leaf_Zn<3000))]
#Mias_a42
#1991.748 Leaf Zn
#174.645 X_Zn 22.815 X_Cd 183.7379 Leaf_Cd 2.062667 CH_Cd 5.581667 CH_Zn
mean(icp$X_Zn[icp$Species=="arenosa"&icp$Population=="Mias"],na.rm=T)
#2089.983
mean(icp$X_Cd[icp$Species=="arenosa"&icp$Population=="Mias"],na.rm=T)
#80.29492
mean(icp$CH_Zn[icp$Species=="arenosa"&icp$Population=="Mias"],na.rm=T)
#239.8335
mean(icp$CH_Cd[icp$Species=="arenosa"&icp$Population=="Mias"],na.rm=T)
#11.16752















####
#Fe#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa"|icp_df$Species=="hybrid",]
icp_df$Species<-droplevels(as.factor(icp_df$Species))
icp_df$Population<-droplevels(as.factor(icp_df$Population))
boxplot(log10(icp_df$Leaf_Fe)~icp_df$Species*icp_df$Population,las=2,ylim=c(1,5),border=c("green","darkblue","darkcyan"),xlab="",ylab=expression(Log[10](Leaf~Fe)*"[ppm]"),cex.lab=1.3,cex.axis=1.2,na.action=na.exclude)

pdf("Leaf_ICP_introgression_boxplots_Fe.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa","hybrid"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df$Leaf_Fe)~icp_df$Species*icp_df$Population,las=2,ylim=c(1,4.5),border=c(rep(c("darkgreen","darkblue","darkcyan"),6)),names=c(rep("",18)),xlab="",ylab=expression(Log[10]*"(Leaf Fe [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:2,4:5,7:11,13:17)],las=2,ylim=c(1,4.5),border=c("darkgreen","darkblue","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10]*"(Leaf Fe [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10(icp_df$Leaf_Fe)~icp_df$Species*icp_df$Population,at=c(1:2,NA,NA,4,NA,5:9,NA,10:12,13:14,NA),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue","darkcyan"),6)))
stripchart(log10(icp_df$Leaf_Fe)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df$Leaf_Fe)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,0.75, 2.4,0.75, xpd = TRUE,lwd=2)
segments(2.6,0.75, 4.4,0.75, xpd = TRUE,lwd=2)
segments(4.6,0.75, 7.4,0.75, xpd = TRUE,lwd=2)
segments(7.6,0.75, 9.4,0.75, xpd = TRUE,lwd=2)
segments(9.6,0.75, 12.4,0.75, xpd = TRUE,lwd=2)
segments(12.6,0.75, 14.4,0.75, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,6,8.5,11,13.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue","darkcyan"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa),italic(Hybrid))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_wohybrids_Fe.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df2$Leaf_Fe)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(1,4.5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Fe [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(1,4.5),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Fe [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10(icp_df2$Leaf_Fe)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Fe)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Fe)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,0.75, 2.4,0.75, xpd = TRUE,lwd=2)
segments(2.6,0.75, 4.4,0.75, xpd = TRUE,lwd=2)
segments(4.6,0.75, 6.4,0.75, xpd = TRUE,lwd=2)
segments(6.6,0.75, 8.4,0.75, xpd = TRUE,lwd=2)
segments(8.6,0.75, 10.4,0.75, xpd = TRUE,lwd=2)
segments(10.6,0.75, 12.4,0.75, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("bottomleft",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

####
#Cu#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa"|icp_df$Species=="hybrid",]
icp_df$Species<-droplevels(as.factor(icp_df$Species))
icp_df$Population<-droplevels(as.factor(icp_df$Population))
boxplot(log10(icp_df$Leaf_Cu)~icp_df$Species*icp_df$Population,las=2,ylim=c(1,5),border=c("green","darkblue","darkcyan"),xlab="",ylab=expression(Log[10](Leaf~Cu)*"[ppm]"),cex.lab=1.3,cex.axis=1.2,na.action=na.exclude)

pdf("Leaf_ICP_introgression_boxplots_Cu.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa","hybrid"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df$Leaf_Cu)~icp_df$Species*icp_df$Population,las=2,ylim=c(0,2),border=c(rep(c("darkgreen","darkblue","darkcyan"),6)),names=c(rep("",18)),xlab="",ylab=expression(Log[10]*"(Leaf Cu [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:2,4:5,7:11,13:17)],las=2,ylim=c(0,2),border=c("darkgreen","darkblue","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10]*"(Leaf Cu [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10(icp_df$Leaf_Cu)~icp_df$Species*icp_df$Population,at=c(1:2,NA,NA,4,NA,5:9,NA,10:12,13:14,NA),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue","darkcyan"),6)))
stripchart(log10(icp_df$Leaf_Cu)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df$Leaf_Cu)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,-0.1, 2.4,-0.1, xpd = TRUE,lwd=2)
segments(2.6,-0.1, 4.4,-0.1, xpd = TRUE,lwd=2)
segments(4.6,-0.1, 7.4,-0.1, xpd = TRUE,lwd=2)
segments(7.6,-0.1, 9.4,-0.1, xpd = TRUE,lwd=2)
segments(9.6,-0.1, 12.4,-0.1, xpd = TRUE,lwd=2)
segments(12.6,-0.1, 14.4,-0.1, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,6,8.5,11,13.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue","darkcyan"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa),italic(Hybrid))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_wohybrids_Cu.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df2$Leaf_Cu)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,2),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Cu [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,2),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Cu [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10(icp_df2$Leaf_Cu)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Cu)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Cu)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,-0.1, 2.4,-0.1, xpd = TRUE,lwd=2)
segments(2.6,-0.1, 4.4,-0.1, xpd = TRUE,lwd=2)
segments(4.6,-0.1, 6.4,-0.1, xpd = TRUE,lwd=2)
segments(6.6,-0.1, 8.4,-0.1, xpd = TRUE,lwd=2)
segments(8.6,-0.1, 10.4,-0.1, xpd = TRUE,lwd=2)
segments(10.6,-0.1, 12.4,-0.1, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()


pdf("Leaf_ICP_introgression_boxplots_wohybrids_Cu_nolog.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(icp_df2$Leaf_Cu~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0,50),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Cu [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0,50),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Cu [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(icp_df2$Leaf_Cu~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(icp_df2$Leaf_Cu[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(icp_df2$Leaf_Cu[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,-3, 2.4,-3, xpd = TRUE,lwd=2)
segments(2.6,-3, 4.4,-3, xpd = TRUE,lwd=2)
segments(4.6,-3, 6.4,-3, xpd = TRUE,lwd=2)
segments(6.6,-3, 8.4,-3, xpd = TRUE,lwd=2)
segments(8.6,-3, 10.4,-3, xpd = TRUE,lwd=2)
segments(10.6,-3, 12.4,-3, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()





####
#Mn#
####
icp_df<-icp[icp$Population=="Buko"|icp$Population=="Kato"|icp$Population=="Mias"|icp$Population=="Piek"|icp$Population=="Zapa"|icp$Population=="Kowa",]
icp_df<-icp_df[icp_df$Species=="halleri"|icp_df$Species=="arenosa"|icp_df$Species=="hybrid",]
icp_df$Species<-droplevels(as.factor(icp_df$Species))
icp_df$Population<-droplevels(as.factor(icp_df$Population))
boxplot(log10(icp_df$Leaf_Mn)~icp_df$Species*icp_df$Population,las=2,ylim=c(1,5),border=c("green","darkblue","darkcyan"),xlab="",ylab=expression(Log[10](Leaf~Mn)*"[ppm]"),cex.lab=1.3,cex.axis=1.2,na.action=na.exclude)

pdf("Leaf_ICP_introgression_boxplots_Mn.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df$Species <- factor(as.factor(icp_df$Species),levels=c("halleri","arenosa","hybrid"))
icp_df$Population <- factor(as.factor(icp_df$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df$Leaf_Mn)~icp_df$Species*icp_df$Population,las=2,ylim=c(0.5,3),border=c(rep(c("darkgreen","darkblue","darkcyan"),6)),names=c(rep("",18)),xlab="",ylab=expression(Log[10]*"(Leaf Mn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:2,4:5,7:11,13:17)],las=2,ylim=c(0.5,3),border=c("darkgreen","darkblue","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue","darkgreen","darkblue","darkcyan","darkgreen","darkblue"),names=c(rep("",14)),xlab="",ylab=expression(Log[10]*"(Leaf Mn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10(icp_df$Leaf_Mn)~icp_df$Species*icp_df$Population,at=c(1:2,NA,NA,4,NA,5:9,NA,10:12,13:14,NA),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue","darkcyan"),6)))
stripchart(log10(icp_df$Leaf_Mn)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df$Leaf_Mn)[icp_df$Species=="halleri"&icp_df$Population=="Kato"&icp_df$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,0.35, 2.4,0.35, xpd = TRUE,lwd=2)
segments(2.6,0.35, 4.4,0.35, xpd = TRUE,lwd=2)
segments(4.6,0.35, 7.4,0.35, xpd = TRUE,lwd=2)
segments(7.6,0.35, 9.4,0.35, xpd = TRUE,lwd=2)
segments(9.6,0.35, 12.4,0.35, xpd = TRUE,lwd=2)
segments(12.6,0.35, 14.4,0.35, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,6,8.5,11,13.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue","darkcyan"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa),italic(Hybrid))),cex=1.5,bty="n",pt.lwd=4)
dev.off()

pdf("Leaf_ICP_introgression_boxplots_wohybrids_Mn.pdf",width=8,height=8,paper="special",pointsize=14)
icp_df2<-icp_df[!icp_df$Species=="hybrid",]
icp_df2<-droplevels(icp_df2)
icp_df2$Species <- factor(as.factor(icp_df2$Species),levels=c("halleri","arenosa"))
icp_df2$Population <- factor(as.factor(icp_df2$Population),levels=c("Buko","Kato","Mias","Piek","Zapa","Kowa"))
par(lwd=1)
par(mar=c(6,5,1,1))

x<-boxplot(log10(icp_df2$Leaf_Mn)~icp_df2$Species*icp_df2$Population,las=2,ylim=c(0.5,3),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Mn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,outwex=1.5)
boxplot(x$stats[,c(1:12)],las=2,ylim=c(0.5,3),border=c("darkgreen","darkblue"),names=c(rep("",12)),xlab="",ylab=expression(Log[10]*"(Leaf Mn [ppm])"),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5)
stripchart(log10(icp_df2$Leaf_Mn)~icp_df2$Species*icp_df2$Population,at=c(1:2,NA,4:12),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("darkgreen","darkblue"),6)))
stripchart(log10(icp_df2$Leaf_Mn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID!="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=17,col=c("darkgreen"))
stripchart(log10(icp_df2$Leaf_Mn)[icp_df2$Species=="halleri"&icp_df2$Population=="Kato"&icp_df2$Sample_ID=="Kato_h09"],at=c(3),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c("darkgreen"))

segments(0.6,0.35, 2.4,0.35, xpd = TRUE,lwd=2)
segments(2.6,0.35, 4.4,0.35, xpd = TRUE,lwd=2)
segments(4.6,0.35, 6.4,0.35, xpd = TRUE,lwd=2)
segments(6.6,0.35, 8.4,0.35, xpd = TRUE,lwd=2)
segments(8.6,0.35, 10.4,0.35, xpd = TRUE,lwd=2)
segments(10.6,0.35, 12.4,0.35, xpd = TRUE,lwd=2)

mtext(c("Buko","Kato","Mias","Piek","Zapa","Kowa"),side=1,line=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),cex=1.5,col=c(rep("red",4),"black","black"))
par(lwd=2)
legend("topright",fill="white",border=c("darkgreen","darkblue"),legend=c(expression(italic(A.~halleri),italic(A.~arenosa))),cex=1.5,bty="n",pt.lwd=4)
dev.off()









