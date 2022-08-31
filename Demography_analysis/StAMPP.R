#Stampp
require(StAMPP)
require(vcfR)
require(adegenet)
# e.g. fst <- stamppFst(geno)

#from Patrick Monnahan: 1.6.21
# convert to genlight 
vcf<-read.vcfR("All5_introgression_4dgsites.LD_Pruned.vcf", verbose = FALSE)
##2.2 Convert vcf into genlight object
#leads to NA for tetraploids
Intro_GL <- vcfR2genlight.tetra(vcf)
locNames(Intro_GL) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(Intro_GL)<-gsub("Zako","Zapa",substr(indNames(Intro_GL),1,4))               # add pop names: here pop names are first 3 chars of ind name

#check
Intro_GL
#186 genotypes,  23,301 binary SNPs, size: 4.5 Mb
# 140445 (3.24 %) missing data

indNames(Intro_GL)
ploidy(Intro_GL)
pop(Intro_GL)

##2.3 Convert genlight to allele frequencies for StAMPP
#Intro_StAMPP<-stamppConvert(Intro_GL, type = "genlight")

############
#   PCA 
############
# run PCA
pca.1 <- glPcaFast(Intro_GL, nf=186) # use the modified function glPcaFast at the end of the file
# proportion of explained variance by first three axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis

colour_pops<-c(rep("red",28)#Outgroup
,rep("green",26)#Buko
,rep("blue",19)#Kato
,rep("orange",19)#Kowa
,rep("purple",34)#Mias
,rep("lightblue",34)#Piek
,rep("yellow",26)#Zapa
)

pdf("PCA_Introgression_woutgroup_LDpruned_MAF_adegenet.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(pca.1$scores[,1],pca.1$scores[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(pca.1$eig[1]/sum(pca.1$eig)*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(pca.1$eig[2]/sum(pca.1$eig)*100,0),"%)"),cex=1.5)
text(pca.1$scores[,1],pca.1$scores[,2],labels=rownames(pca.1$scores), cex= 0.7)
dev.off()

pdf("PCA_Introgression_woutgroup_LDpruned_MAF_adegenet_PC1_3.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,8))
par(mgp=c(3.5,0.75,0))
plot(pca.1$scores[,1],pca.1$scores[,3], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("bottomright",fill=unique(colour_pops),legend=c("Outgroup","Buko","Kato","Kowa","Mias","Piek","Zapa"),cex=1.3,border=c("black"),inset=c(-0.45,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(pca.1$eig[1]/sum(pca.1$eig)*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC3 (",round(pca.1$eig[3]/sum(pca.1$eig)*100,0),"%)"),cex=1.5)
text(pca.1$scores[,1],pca.1$scores[,3],labels=rownames(pca.1$scores), cex= 0.7)
dev.off()

############################
#     K-means clustering
############################
## find clusters = K-means clustering
grp <- find.clusters(Intro_GL, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000)
#186 PCs retained
#7 clusters

# save the grouping
write.table(grp$grp, file="K_means_introgression.txt", sep="\t", quote=F, col.names=F)

############################
#     K-means clustering no outgroup
############################
## find clusters = K-means clustering
pca.2 <- glPcaFast(Intro_GL[29:186], nf=158) # use the modified function glPcaFast at the end of the file
grp <- find.clusters(Intro_GL[29:186], max.n.clust=20, glPca = pca.2, perc.pca = 100, n.iter=1e6, n.start=1000)
#158 PCs retained
#2 clusters

# save the grouping
write.table(grp$grp, file="K_means_introgression_nooutgroup.txt", sep="\t", quote=F, col.names=F)

###################################################
####   distance-based analyses  #######
###################################################

### Calculate Nei's distances between individuals/pops
# ---------------------------------------------------
aa.D.ind <- stamppNeisD(Intro_GL, pop = FALSE)  # Nei's 1972 distance between indivs
stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance.phy.dst") # export matrix - for SplitsTree
aa.D.pop <- stamppNeisD(Intro_GL, pop = TRUE)   # Nei's 1972 distance between pops
stamppPhylip(aa.D.pop, file="aa.pops_Neis_distance.phy.dst") # export matrix - for SplitsTree

aa.D.ind <- stamppNeisD(Intro_GL[29:186], pop = FALSE)  # Nei's 1972 distance between indivs
stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance.phy.MAF.dst") # export matrix - for SplitsTree

### create the dist objects
colnames(aa.D.ind) <- rownames(aa.D.ind) 
aa.D.ind.dist<-as.dist(aa.D.ind, diag=T)
attr(aa.D.ind.dist, "Labels")<-rownames(aa.D.ind)          # name the rows of a matrix  

colnames(aa.D.pop) <- rownames(aa.D.pop) 
aa.D.pop.dist<-as.dist(aa.D.pop, diag=T)
attr(aa.D.pop.dist, "Labels")<-rownames(aa.D.pop)          # name the rows of a matrix  

require(openxlsx)
write.xlsx(aa.D.ind,"NeisD_introgression_ind.xlsx",row.names=T)
write.xlsx(aa.D.pop,"NeisD_introgression_pop.xlsx",row.names=T)

require(dartR)

Intro_GI<-gl2gi(Intro_GL)
save(Intro_GI, file="Intro_GI.rdata")









##################### 
# MODIFIED FUNCTIONS

# a function for conversion from vcfR object to genlight in tetraploids
vcfR2genlight.tetra <- function (x, n.cores = 1) 
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

## -------------------------------------------------------
### a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}



Intro_het<-read.table("Heterozygosity_introgression_stats5",header=F)
colnames(Intro_het)<-c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions","nIndels","average depth","nSingletons","nHapRef","nHapAlt","nMissing")
percHet=Intro_het$nHets/(Intro_het$nNonRefHom+Intro_het$nHets)*100
Intro_het$percHet<-percHet
require(openxlsx)
write.xlsx(Intro_het,"Heterozygosity_Introgression_bcftools.xlsx",sheetName="Introgression",overwrite=T)




#Nei's D per population and ploidy

pop(Intro_GL)<-paste(gsub("Zako","Zapa",substr(indNames(Intro_GL),1,4)),ploidy(Intro_GL),sep="_")

aa.D.pop <- stamppNeisD(Intro_GL, pop = TRUE)   # Nei's 1972 distance between pops
stamppPhylip(aa.D.pop, file="aa.pops_Neis_distance.phy.dst") # export matrix - for SplitsTree

colnames(aa.D.pop) <- rownames(aa.D.pop) 
aa.D.pop.dist<-as.dist(aa.D.pop, diag=T)
attr(aa.D.pop.dist, "Labels")<-rownames(aa.D.pop)          # name the rows of a matrix  

require(openxlsx)
write.xlsx(as.data.frame(aa.D.pop),"NeisD_introgression_pop_ploidy.xlsx",row.names=T,overwrite=T)


