require(polyRAD)
require(VariantAnnotation)
mydata<-VCF2RADdata("All5_introgression_4dgsites.LD_Pruned.vcf",possiblePloidies=list(2,4),expectedLoci=23301,expectedAlleles=46602,phaseSNPs=FALSE,
min.ind.with.reads=1,min.ind.with.minor.allele=2,contamRate=0.001,samples = VariantAnnotation::samples(VariantAnnotation::scanVcfHeader("All5_introgression_4dgsites.LD_Pruned.vcf")))
#23188 loci
#Next we can check for markers that are behaving in a non-Mendelian fashion. If we are expecting diploid segregation, all markers should show a Hind/HE value of 0.5 or less. (For an autopolyploid, the expected value is (ploidy-1)/ploidy.)

myhindhe <- HindHe(mydata)
myhindheByLoc <- colMeans(myhindhe, na.rm = TRUE)
hist(myhindheByLoc, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus")
abline(v = 0.5, col = "blue", lwd = 2)

#Hind/HE for filtering markers and individuals

#Some markers may behave in a non-Mendelian fashion due to misalignments, amplification bias, presence-absence variation, or other issues. In addition to filtering out problematic markers, you may also want to confirm that all individuals in the dataset are well-behaved.

#The Hind/HE statistic (Clark et al. 2020) helps to filter such markers and individuals. In a mapping population it can be run using the HindHeMapping function, which requires a single ploidy to be input, along with the mapping population design. In a natural population or diversity panel, the HindHe function can be used. HindHe should also be used for mapping populations in which the most recent generation was created by random intermating among all progeny. In all cases, I recommend running HindHe or HindHeMapping before running TestOverdispersion or any of the genotype calling functions, as demonstrated in the previous sections.

TotDepthT <- rowSums(mydata$locDepth)

#TotDepthT is a vector showing the total read depth at each locus.
#To investigate individuals, we can take the row means of the matrix:

myHindHeByInd <- rowMeans(myhindhe, na.rm = TRUE)
hist(myHindHeByInd, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by Ind")
abline(v = 0.5, col = "blue", lwd = 2)

pdf("Hind_He_subset_min2indminorallele_All5.pdf",width=20,height=20,paper="special",pointsize=10)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
HindHeplot<-plot(myHindHeByInd)
text(158,myHindHeByInd[158],names(myHindHeByInd)[158])
text(97,myHindHeByInd[97],names(myHindHeByInd)[97])
text(196,myHindHeByInd[196],names(myHindHeByInd)[196])
text(67,myHindHeByInd[67],names(myHindHeByInd)[67])
text(99,myHindHeByInd[99],names(myHindHeByInd)[99])
text(52,myHindHeByInd[52],names(myHindHeByInd)[52])
text(50,myHindHeByInd[196],names(myHindHeByInd)[50])
dev.off()

pdf("Hind_He_min2indminorallele_All5.pdf",width=20,height=20,paper="special",pointsize=10)
par(mar=c(5,5,1, 8))
par(mgp=c(3.5,0.75,0))
HindHeplot<-plot(myHindHeByInd)
text(1:205,myHindHeByInd,names(myHindHeByInd))
dev.off()

Genotyping_ploidy4<-read.table("../Total_samples_introgression_4.list",sep=",",header=F)
Genotyping_ploidy2<-read.table("../Total_samples_introgression_2.list",sep=",",header=F)
GT_ploidy4<-cbind(Genotyping_ploidy4,rep(4,length(Genotyping_ploidy4)))
colnames(GT_ploidy4)<-c("Sample","Ploidy")
GT_ploidy2<-cbind(Genotyping_ploidy2,rep(2,length(Genotyping_ploidy2)))
colnames(GT_ploidy2)<-c("Sample","Ploidy")
Missing<-data.frame(c("Kowa_h07","Kowa_h41","Kowa_h42","Kowa_h43","Piek_a20_1","Piek_a20_2","Piek_h14","Piek_h6"),c(2,2,2,2,4,4,4,4))
colnames(Missing)<-c("Sample","Ploidy")


Genotyping_ploidy<-rbind(GT_ploidy4,GT_ploidy2,Missing)
Ploidy<-merge(as.data.frame(names(myHindHeByInd)),Genotyping_ploidy,by.x="names(myHindHeByInd)",by.y="Sample",all.x=T)
colnames(Ploidy)<-c("Sample","Ploidy")
colour_pops<-c(rep("blue",28)#Outgroup
,rep("orange",26)#Buko
,rep("magenta",19)#Kato
,rep("black",19)#Kowa
,rep("red",34)#Mias
,rep("purple",34)#Piek
,rep("grey",26)#Zapa
)


Ploidy$Ploidy[1:28]<-2
Ploidy$pop<-substr(Ploidy$Sample,1,4)
Ploidy$pop[Ploidy$pop=="Zako"]<-"Zapa"


pdf("Hind_He_All5_paper.pdf",width=20,height=20,paper="special",pointsize=15)
par(mar=c(5,5,1, 8))
par(mgp=c(2.5,0.75,0))
HindHeplot<-plot(myHindHeByInd,pch=c(17,19)[as.factor(Ploidy$Ploidy[match(names(myHindHeByInd),Ploidy$Sample)])],las=1,xlab="Plant number",cex=1.5)
dev.off()

myHindHeByInd2<-myHindHeByInd[-c(1:28)]
myHindHeByInd2<-myHindHeByInd2[order(Ploidy$pop[match(names(myHindHeByInd2),Ploidy$Sample)],Ploidy$Ploidy[match(names(myHindHeByInd2),Ploidy$Sample)])]
myHindHeByInd2<-myHindHeByInd2[c(46:64,133:158,27:45,1:26,99:132,65:98)]
colour_pops2<-c("orange","magenta","black","red","purple","grey")[as.factor(Ploidy$pop[match(names(myHindHeByInd2),Ploidy$Sample)])]

pdf("Hind_He_All5_paper_wooutgroup.pdf",width=12,height=12,paper="special",pointsize=25)
par(lwd=2)
par(mar=c(5,5,4,1))
par(mgp=c(2.5,0.75,0))
HindHeplot<-plot(myHindHeByInd2,pch=c(17,19)[as.factor(Ploidy$Ploidy[match(names(myHindHeByInd2),Ploidy$Sample)])],col=colour_pops2,ylab=expression(H[ind]/H[E]),cex.lab=1.5,cex.axis=1.3,las=1,xlab="Plant identity (no)")
#text(40:50,myHindHeByInd2[40:50],names(myHindHeByInd2[40:50]))
dev.off()


#expected value under Hardy-Weinberg Equilibrium. This is ploidy-1/ploidy, e.g. 0.5 for diploids and 0.75 for tetraploids. Since there is some population structure, most individuals show a lower value. However, some interspecific hybrids have values higher than expected. We can also see that it is fairly easy to distinguish diploids and tetraploids. This method is not a replacement for flow cytometry, but can complement it if some minority of samples in the dataset have unknown ploidy.

require(openxlsx)
write.xlsx(as.data.frame(myHindHeByInd),"HindHeByInd_Introgression_min2indminorallele_All5.xlsx",row.names=T)
