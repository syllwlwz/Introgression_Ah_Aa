require(openxlsx)
data<-read.xlsx("F:/Introgression/ICP_arenosa2/Vero_leaf_ICP/ICP_tolexp_Veronica_30102017_LS.xlsx",1)
summary.factor(data$Group)
#   A     
# 80      
icp<-data[!is.na(data$Group),]
icp$Species <- relevel(as.factor(icp$Species), "halleri")
icp$SoilType <- relevel(as.factor(icp$SoilType), "Mias")


pdf("ICP_Vero_tolexp_paper_Zn.pdf",width=10,height=12,paper="special",pointsize=28)
par(lwd=2.5)
par(mar=c(5,7,5,1))
par(mgp=c(4.5,0.75,0))
x<-boxplot(log10(icp$Zn)~icp$Pop*icp$Species*icp$SoilType,las=2,ylim=c(1,5),border=c("red","black"),names=c(rep("",8)),xlab="",ylab=expression(Leaf~Zn~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,plot=F)
boxplot(x$stats[,c(1:8)],las=2,ylim=c(1,5),border=c("red","black"),names=c(rep("",8)),xlab="",ylab=expression(Leaf~Zn~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n")
abline(h=log10(3000),col="grey40",lwd=2.5)
stripchart(log10(icp$Zn)~icp$Pop*icp$Species*icp$SoilType,at=c(1:8),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("red","black"),4)))
axis(side=2,at=c(1:5),labels=c("10","100","1,000","10,000","100,000"),las=1,cex.axis=1.3)
segments(0.6,0.75, 2.4,0.75, xpd = TRUE)
segments(2.6,0.75, 4.4,0.75, xpd = TRUE)
segments(4.6,0.75, 6.4,0.75, xpd = TRUE)
segments(6.6,0.75, 8.4,0.75, xpd = TRUE)
mtext(c(expression(bold(italic(Ah)),bold(italic(Aa)))),side=1,line=1,at=c(1.5,3.5,5.5,7.5),cex=1.5)
segments(0.6,0.25, 4.4,0.25, xpd = TRUE)
segments(4.6,0.25, 8.4,0.25, xpd = TRUE)
mtext(c("Mias soil","NM soil"),side=1,line=3,at=c(2.5,6.5),cex=1.5,font=2)
par(lwd=2)
legend(-4,6.2,fill="white",border=c("red","black"),legend=c("Mias origin","Zapa origin"),cex=1.5,bty="n",title="Plants",horiz=T,xpd=T)
dev.off()


pdf("ICP_Vero_tolexp_paper_Cd.pdf",width=10,height=10,paper="special",pointsize=28)
par(lwd=2.5)
par(mar=c(5,7,1,1))
par(mgp=c(3.5,0.75,0))
x<-boxplot(log10(icp$Cd+1)~icp$Pop*icp$Species*icp$SoilType,las=2,ylim=c(0,3),border=c("red","black"),names=c(rep("",8)),xlab="",ylab=expression(Leaf~Cd~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,plot=F)
boxplot(x$stats[,c(1:3,NA,5:8)],at=c(1:8),las=2,ylim=c(0,3),border=c("red","black"),names=c(rep("",8)),xlab="",ylab=expression(Leaf~Cd~conc.~(µg~g^-1~DW)),xaxt="n",cex.lab=1.5,cex.axis=1.3,col="white",outwex=1.5,yaxt="n")
abline(h=log10(100),col="grey40",lwd=2.5)
stripchart(log10(icp$Cd+1)~icp$Pop*icp$Species*icp$SoilType,at=c(1:3,NA,5:8),vertical=TRUE,method="jitter",add=TRUE,pch=20,col=c(rep(c("red","black"),4)))
axis(side=2,at=c(1:4),labels=c("10","100","1,000","10,000"),las=1,cex.axis=1.3)
segments(0.6,-0.2, 2.4,-0.2, xpd = TRUE)
segments(2.6,-0.2, 4.4,-0.2, xpd = TRUE)
segments(4.6,-0.2, 6.4,-0.2, xpd = TRUE)
segments(6.6,-0.2, 8.4,-0.2, xpd = TRUE)
mtext(c(expression(bold(italic(Ah)),bold(italic(Aa)))),side=1,line=1,at=c(1.5,3.5,5.5,7.5),cex=1.5)
segments(0.6,-0.6, 4.4,-0.6, xpd = TRUE)
segments(4.6,-0.6, 8.4,-0.6, xpd = TRUE)
mtext(c("Mias soil","NM soil"),side=1,line=3,at=c(2.5,6.5),cex=1.5,font=2)
par(lwd=2)
#legend("topright",fill="white",border=c("red","black"),legend=c("Mias origin","Zapa origin"),cex=1.5,bty="n",title="Plants")
dev.off()














############
#Statistics#
############
data<-read.xlsx("F:/Introgression/ICP_arenosa2/Vero_leaf_ICP/ICP_tolexp_Veronica_30102017_LS.xlsx",1)
summary.factor(data$Group)
#   A     
# 80        
icp<-data[!(is.na(data$Group)|is.na(data$Zn)),]
icp$Species <- relevel(as.factor(icp$Species), "halleri")
icp$SoilType <- relevel(as.factor(icp$SoilType), "Mias")
icp$Pop <- as.factor(icp$Pop)
icp<-droplevels(icp)

require (car)
require (MASS)
require (vegan)
require(multcompView)

#...............
#Equal variances
#...............

leveneTest(log10(icp$Zn)~icp$Pop*icp$Species*icp$SoilType)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value    Pr(>F)    
#group  7  0.9457  0.468
#      77 
    
#Levene test: If p<0.05, the results of ANOVA are less reliable. There is no equivalent test but comparing the 
#p-values from the ANOVA with 0.01 instead of 0.05 is acceptable.
#Levene is less strict because relies less on normal distrib. of variables

###############
#ANOVA
###############

anova<-aov(log10(icp$Zn)~icp$Pop*icp$Species*icp$SoilType,data=icp)
hist(anova$residuals)
shapiro.test(anova$residuals) 
#W = 0.7642, p-value = 5.297e-10
#can not use Anova

############################
#T-test#
############################
require (car)
require (MASS)
require (vegan)
require(multcompView)
require(PairedData)

data<-icp
data$Comb<-paste(icp$Pop, icp$Species, icp$SoilType,sep=":")
data$Comb<-as.factor(data$Comb)
data$Comb<-droplevels(data$Comb)
Test_df<-data.frame(matrix(ncol=3,nrow=0))
colnames(Test_df)<-c("Test_combination","p-value","Statistic")
for (i in levels(data$Comb))
	{for (j in levels(data$Comb))
		{if (i==j)
			{next}
		else if (shapiro.test(data$Zn[data$Comb==i|data$Comb==j])$p.value<=0.05)
			{if (leveneTest(data$Zn[data$Comb==i|data$Comb==j]~data$Comb[data$Comb==i|data$Comb==j])[1,3]<=0.05)
				{Comb<-paste(i,j,sep="-")
				Test<-yuen.t.test(data$Zn[data$Comb==i],data$Zn[data$Comb==j])$p.value
				Statistic<-"Yuen"
				Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
				}
			else
				{Comb<-paste(i,j,sep="-")
				Test<-wilcox.test(data$Zn[data$Comb==i],data$Zn[data$Comb==j],exact=T)$p.value
				Statistic<-"Wilcox"
				Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
				}
			next
			}
		else if (leveneTest(data$Zn[data$Comb==i|data$Comb==j]~data$Comb[data$Comb==i|data$Comb==j])[1,3]<=0.05)
			{Comb<-paste(i,j,sep="-")
			Test<-t.test(data$Zn[data$Comb==i],data$Zn[data$Comb==j])$p.value
			Statistic<-"Welch_t_test"
			Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
			next
			}
		else
			{Comb<-paste(i,j,sep="-")
			Test<-t.test(data$Zn[data$Comb==i],data$Zn[data$Comb==j],var.equal=T)$p.value
			Statistic<-"T_test"
			Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
			}
		}
	}

#Warnmeldung:
#In wilcox.test.default(data$Root_biomass[data$Comb == i], data$Root_biomass[data$Comb ==  :
#  kann bei Bindungen keinen exakten p-Wert Berechnen
#wenn exakt gleiche Werte vorkommen

Test_df[,2]<-unname(Test_df[,2],force=TRUE)
require(qvalue)
Test_df$p_adj<-qvalue(as.numeric(as.character(Test_df[,2])),lambda=0)$qvalue
colnames(Test_df)<-c("Test_combination","p-value","Statistic","p adj")
rownames(Test_df)<-Test_df[,1]

multcompLetters(extract_p(Test_df))
multcompLetters2(Zn~Comb,extract_p(Test_df),data=data)
#   Zapa:halleri:Mias    Mias:halleri:Mias    Mias:arenosa:Mias 
#                 "a"                  "a"                  "b" 
#Zapa:halleri:Control Mias:halleri:Control Mias:arenosa:Control 
#                 "c"                  "c"                  "d" 
#Zapa:arenosa:Control 
#                 "e"

write.xlsx(Test_df,"T_tests_tolexp_Zn.xlsx",sheetName="Zn",row.names=F)
#a * means data was neither normally distributed nor had equal variances. Welch t-test was employed but treat with caution



data<-icp
data$Comb<-paste(icp$Pop, icp$Species, icp$SoilType,sep=":")
data$Comb<-as.factor(data$Comb)
data$Comb<-droplevels(data$Comb)
Test_df<-data.frame(matrix(ncol=3,nrow=0))
colnames(Test_df)<-c("Test_combination","p-value","Statistic")
for (i in levels(data$Comb))
	{for (j in levels(data$Comb))
		{if (i==j)
			{next}
		else if (shapiro.test(data$Cd[data$Comb==i|data$Comb==j])$p.value<=0.05)
			{if (leveneTest(data$Cd[data$Comb==i|data$Comb==j]~data$Comb[data$Comb==i|data$Comb==j])[1,3]<=0.05)
				{Comb<-paste(i,j,sep="-")
				Test<-yuen.t.test(data$Cd[data$Comb==i],data$Cd[data$Comb==j])$p.value
				Statistic<-"Yuen"
				Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
				}
			else
				{Comb<-paste(i,j,sep="-")
				Test<-wilcox.test(data$Cd[data$Comb==i],data$Cd[data$Comb==j],exact=T)$p.value
				Statistic<-"Wilcox"
				Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
				}
			next
			}
		else if (leveneTest(data$Cd[data$Comb==i|data$Comb==j]~data$Comb[data$Comb==i|data$Comb==j])[1,3]<=0.05)
			{Comb<-paste(i,j,sep="-")
			Test<-t.test(data$Cd[data$Comb==i],data$Cd[data$Comb==j])$p.value
			Statistic<-"Welch_t_test"
			Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
			next
			}
		else
			{Comb<-paste(i,j,sep="-")
			Test<-t.test(data$Cd[data$Comb==i],data$Cd[data$Comb==j],var.equal=T)$p.value
			Statistic<-"T_test"
			Test_df<-rbind(Test_df,cbind(Comb,Test,Statistic))
			}
		}
	}

#Warnmeldung:
#In wilcox.test.default(data$Root_biomass[data$Comb == i], data$Root_biomass[data$Comb ==  :
#  kann bei Bindungen keinen exakten p-Wert Berechnen
#wenn exakt gleiche Werte vorkommen

Test_df[,2]<-unname(Test_df[,2],force=TRUE)
require(qvalue)
Test_df$p_adj<-qvalue(as.numeric(as.character(Test_df[,2])),lambda=0)$qvalue
colnames(Test_df)<-c("Test_combination","p-value","Statistic","p adj")
rownames(Test_df)<-Test_df[,1]

multcompLetters(extract_p(Test_df))
multcompLetters2(Cd~Comb,extract_p(Test_df),data=data)
#   Zapa:halleri:Mias    Mias:halleri:Mias    Mias:arenosa:Mias 
#                 "a"                  "a"                  "b" 
#Zapa:halleri:Control Mias:halleri:Control Zapa:arenosa:Control 
#                 "c"                  "d"                  "e" 
#Mias:arenosa:Control 
#                 "e"

write.xlsx(Test_df,"T_tests_tolexp_Cd.xlsx",sheetName="Cd",row.names=F)
#a * means data was neither normally distributed nor had equal variances. Welch t-test was employed but treat with caution


