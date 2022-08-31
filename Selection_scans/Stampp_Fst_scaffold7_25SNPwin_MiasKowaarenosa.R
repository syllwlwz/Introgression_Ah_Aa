#Stampp
options(java.parameters = "-Xmx12g")

require(StAMPP)
require(vcfR)
require(adegenet)

#from Patrick Monnahan: 1.6.21
##################### 
# MODIFIED FUNCTIONS

# a function for conversion from vcfR object to genlight in tetraploids
vcfR2genlight.tetra <- function (x, n.cores = 4) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}

# convert to genlight 
#vcf<-read.vcfR("MiasKowaarenosa.vcf", verbose = FALSE)
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_1",],file="MiasKowaarenosa_scaffold1.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_2",],file="MiasKowaarenosa_scaffold2.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_3",],file="MiasKowaarenosa_scaffold3.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_4",],file="MiasKowaarenosa_scaffold4.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_5",],file="MiasKowaarenosa_scaffold5.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_6",],file="MiasKowaarenosa_scaffold6.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_7",],file="MiasKowaarenosa_scaffold7.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_8",],file="MiasKowaarenosa_scaffold8.vcf")
#write.vcf(vcf[vcf@ fix[,1]=="scaffold_9",],file="MiasKowaarenosa_scaffold9.vcf")

vcf<-read.vcfR("MiasKowaarenosa_scaffold7.vcf", verbose = FALSE)

##2.2 Convert vcf into genlight object
Intro_GL <- vcfR2genlight.tetra(vcf)
locNames(Intro_GL) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(Intro_GL)<-substr(indNames(Intro_GL),1,4)               # add pop names: here pop names are first 3 chars of ind name
pop(Intro_GL)<-paste(gsub("Kowa","Kowa",substr(indNames(Intro_GL),1,4)),ploidy(Intro_GL),sep="_")

#check
Intro_GL
#69 genotypes,  1,102,246 binary SNPs
#4859345 (6.39 %) missing data

indNames(Intro_GL)
ploidy(Intro_GL)
pop(Intro_GL)

##2.3 Convert genlight to allele frequencies for StAMPP

data.freq<-stamppConvert(Intro_GL,type = "genlight")

Fst <- stamppFst(data.freq, nclusters = 2)
Fst$Fsts
#           Mias_4    Mias_2    Kowa_2 Kowa_4
#Mias_4         NA        NA        NA     NA
#Mias_2 0.46938788        NA        NA     NA
#Kowa_2 0.45632526 0.1295224        NA     NA
#Kowa_4 0.07857803 0.5633260 0.5464411     NA

#           Mias_4    Mias_2    Kowa_2 Kowa_4
#Mias_4         NA        NA        NA     NA
#Mias_2 0.46938788        NA        NA     NA
#Kowa_2 0.45632526 0.1295224        NA     NA
#Kowa_4 0.07857803 0.5633260 0.5464411     NA

#subset
MiasKowaarenosa1<-data.freq[data.freq$pop.names=="Mias_4"|data.freq$pop.names=="Kowa_4",]

NoNAs<-colSums(is.na(MiasKowaarenosa1[,6:ncol(MiasKowaarenosa1)]))
test2<-colSums(MiasKowaarenosa1[,6:ncol(MiasKowaarenosa1)],na.rm=T)

MiasKowaarenosa2<-MiasKowaarenosa1[,6:ncol(MiasKowaarenosa1)]
MiasKowaarenosa3<-MiasKowaarenosa2[,(test2!=0)&(nrow(MiasKowaarenosa2)>(test2+NoNAs))]
#29 obs. of  865219 variables, before 1102246 variables
MiasKowaarenosa<-cbind(MiasKowaarenosa1[,1:5],MiasKowaarenosa3)

Fst_win<-matrix()
count=1
for (i in seq(6,ncol(MiasKowaarenosa)-24,25))
	{Fst_win[count] <- stamppFst(MiasKowaarenosa[,c(1:5,i:(i+24))], nclusters = 2)$Fsts[2,1]   # Nei's 1972 distance between pops
	count=count+1
	}


Fst_df<-data.frame(colnames(MiasKowaarenosa)[seq(6,ncol(MiasKowaarenosa)-24,25)],colnames(MiasKowaarenosa)[seq(30,ncol(MiasKowaarenosa),25)],Fst_win)
require(tidyr)
Fst_df2<-separate(Fst_df,1,sep=11,into=c("Scaffold","Start_position"))
Fst_df3<-separate(Fst_df2,3,sep=11,into=c("Scaffold2","End_position"))
Fst_df4<-Fst_df3[,c(1:2,4:5)]
Fst_df4$Scaffold<-"Scaffold_7"
write.table(Fst_df4,"Fst_Mias_Kowa_arenosa_scaffold7_25SNPwins.table",sep="\t",row.names=F,quote=F)

