######################################################################
#Analysis with potential promoter region before gene of 2 kb included#
######################################################################

#################################################
#load metrics and assign headers where necessary#
#################################################
filelist = list.files(pattern = glob2rx("Fst_Mias_Kowa_arenosa_scaffold*_25SNPwins.table")) 
require(gtools)
filelist <- mixedsort(filelist)
fst_25 = do.call(rbind,lapply(filelist,function(fn)read.delim(fn,header=T, sep="\t"))) #merge all count files

#################################
#kill windows bigger than 100 kb#
#################################
fst_25<-fst_25[fst_25$End_position-fst_25$Start_position<100000,]
#31 windows excluded
#248441 windows left

#############################################
#calculate average and median of window size#
#############################################
mean(fst_25$End_pos-fst_25$Start_pos)
#671.7148
median(fst_25$End_pos-fst_25$Start_pos)
#314

#quantiles
quantile(fst_25$End_pos-fst_25$Start_pos,probs=c(0.01,0.1,0.25,0.75,0.9,0.99))
#     1%     10%     25%     75%     90%     99% 
#122  183  234  453  838 7998

#########################################
#find overlaps between windows and genes#
#########################################

library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

## Load GFF files and add promotor region to genes

lyr_txdm = makeTxDbFromGFF(file ="D:/Lyrata/Alyrata_384_v2.1.gene.gff3", format = c("gff3"), organism = "Arabidopsis lyrata")

#lyr_txdm
lyrGenes=genes(lyr_txdm)
names(lyrGenes)<-NULL

#Define a function that allows the addition of promotor regions
extend <- function(x, upstream=0, downstream=0) #from https://support.bioconductor.org/p/78652/
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

# Use the function to extend promotor region 2kb
  # Upstream=2000 adds 2kb upstream from transcription start site
  # Downstream=0 adds 0bp downstream from the gene.
LyrGenes_plusPromot = extend(lyrGenes,upstream=2000, downstream=0)

options(java.parameters = "-Xmx12000m")
require("openxlsx")

Thal_description<-read.xlsx("E:/Lyrata_RBH/Lyr_TAIR_Mapman_descript_2021_RBH_OBH.xlsx")

UtesMetal=read.xlsx("E:/Introgression/metal_homeostasis_2020_v2.xlsx",2)
UtesMetal$AGI_Number<-toupper(as.character(UtesMetal$AGI_Number))
UtesMetal<-UtesMetal[,c(2:4,6,8,9,11,14:16)]

###########################
#select 1% outlier windows#
###########################
Perc1Fst25<-fst_25[fst_25$Fst>=quantile(fst_25$Fst_win,0.99)[[1]],]

### Convert list of significant markers to a "Genomic Ranges" table. 
Perc1Fst25_GRange<-GRanges(seqnames=tolower(Perc1Fst25$Scaffold),ranges=IRanges(start=Perc1Fst25$Start_position,end=Perc1Fst25$End_position))

values(Perc1Fst25_GRange)<-Perc1Fst25[,4]

## Merge Significant Regions with Gene list in the halleri Genome
Perc1Fst25_hallerigenes=mergeByOverlaps(Perc1Fst25_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1Fst25_hallerigenes_df=data.frame(as.data.frame(Perc1Fst25_hallerigenes$Perc1Fst25_GRange),as.data.frame(Perc1Fst25_hallerigenes$LyrGenes_plusPromot))
colnames(Perc1Fst25_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Fst25_hallerigenes_df$gene_start<-ifelse(Perc1Fst25_hallerigenes_df$gene_strand=="+",Perc1Fst25_hallerigenes_df$gene_start+2000,Perc1Fst25_hallerigenes_df$gene_start)
Perc1Fst25_hallerigenes_df$gene_end<-ifelse(Perc1Fst25_hallerigenes_df$gene_strand=="-",Perc1Fst25_hallerigenes_df$gene_end-2000,Perc1Fst25_hallerigenes_df$gene_end)
Perc1Fst25_hallerigenes_df$gene_size<-Perc1Fst25_hallerigenes_df$gene_size-2000

######                                                                     #####
###### REGIONS NOT OVERLAPPING GENES WILL BE LOST FROM LIST AT THIS POINT  #####
######                                                                     #####

Perc1Fst25_hallerigenes_df_w_orthogroup= merge(Thal_description,Perc1Fst25_hallerigenes_df,by.y="Gene",by.x="Alyr_ID", all.y=TRUE)
Perc1Fst25_hallerigenes_df_w_MM_lists<-merge(Perc1Fst25_hallerigenes_df_w_orthogroup,UtesMetal,all.x=T,by.x="Ath_ID",by.y="AGI_Number")

options(java.parameters = "-Xmx12000m")

require(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Fst_25SNPwin_MiasKowa")

writeData(wb, "Fst_25SNPwin_MiasKowa", Perc1Fst25_hallerigenes_df_w_MM_lists, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Selected_genes_Fst_introgression_1_MiasKowa.xlsx", overwrite = TRUE)

nrow(Perc1Fst25_hallerigenes_df_w_MM_lists)
#2767