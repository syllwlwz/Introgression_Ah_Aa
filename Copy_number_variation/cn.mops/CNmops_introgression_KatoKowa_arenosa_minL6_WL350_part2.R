#Cnmops
options(java.parameters = "-Xmx12000m")

require("openxlsx")

KaKoar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KaKoar.xlsx",1)
KKar<-read.xlsx("GenesOverlappingSignificantMarkers_ThalianaOrthologs_MapManDescriptions_MetalOverlaps_KoKaar.xlsx",1)

KaKoar_unique<-KaKoar[!duplicated(KaKoar$Lyr_Gene,KaKoar$Copy_start,KaKoar$Copy_end),]
KKar_unique<-KKar[!duplicated(KKar$Lyr_Gene,KKar$Copy_start,KKar$Copy_end),]

KaKoar_unique$Lyr_Gene<-as.character(KaKoar_unique$Lyr_Gene)
KKar_unique$Lyr_Gene<-as.character(KKar_unique$Lyr_Gene)
KaKoar_unique$CN_class<-as.character(KaKoar_unique$CN_class)
KKar_unique$CN_class<-as.character(KKar_unique$CN_class)

KaKoaronly<-KaKoar_unique[!((KaKoar_unique$Lyr_Gene%in%KKar_unique$Lyr_Gene)&(KaKoar_unique$CN_class==KKar_unique$CN_class[match(KaKoar_unique$Lyr_Gene,KKar_unique$Lyr_Gene)])),]
KKaronly<-KKar_unique[!((KKar_unique$Lyr_Gene%in%KaKoar_unique$Lyr_Gene)&(KKar_unique$CN_class==KaKoar_unique$CN_class[match(KKar_unique$Lyr_Gene,KaKoar_unique$Lyr_Gene)])),]

write.table(KaKoaronly,"KaKoar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")
write.table(KKaronly,"KoKaar_cnvs_only_strict_minL6_WL350.table",row.names=F,sep="\t")











