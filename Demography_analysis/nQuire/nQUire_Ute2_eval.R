nQuire_Ute2<-read.table("nQuire_lrdmodels_introgression_Ute2.txt",sep="\t",header=F)
colnames(nQuire_Ute2)<-c("file","free","dip","tri","tet","d_dip","d_tri","d_tet")
nQuire_Ute2$file<-gsub("/prj/pflaphy-gscan/nQuire_Ute2/","",nQuire_Ute2$file)

require(openxlsx)
flow_cytometry<-read.xlsx("../../Sequ_flow_cytometry_Vero.xlsx",1)
require(stringr)

test<-ifelse(str_count(flow_cytometry$Sample_ID,pattern="_")>1,gsub("_","",flow_cytometry$Sample_ID),flow_cytometry$Sample_ID)
test2<-gsub("Zapa3","Zapa",test)
test3<-gsub("Zapa4","Zapa",test2)
flow_cytometry$Sample_ID<-test3

nQuire_Ute2$Estimated_ploidy<-colnames(nQuire_Ute2[,c(6:8)])[apply(nQuire_Ute2[,c(6:8)],1,which.min)]
#nQuire_depth$Estimated_ploidy[!is.na(nQuire_depth$d_tet)]<-colnames(nQuire_depth[!is.na(nQuire_depth$d_tet),c(6:8)])[apply(nQuire_depth[!is.na(nQuire_depth$d_tet),c(6:8)],1,which.min)]
#nQuire_depth$Estimated_ploidy[is.na(nQuire_depth$d_tet)]<-NA

nQuire_Ute2$Estimated_ploidy[nQuire_Ute2$Estimated_ploidy=="d_dip"]<-2
nQuire_Ute2$Estimated_ploidy[nQuire_Ute2$Estimated_ploidy=="d_tri"]<-3
nQuire_Ute2$Estimated_ploidy[nQuire_Ute2$Estimated_ploidy=="d_tet"]<-4

#nQuire_Ute2$file<-as.character(nQuire_Ute2$file)
#test<-ifelse((str_count(nQuire_Ute2$file,pattern="_")>1)&grep('_[a-z]$',nQuire_Ute2$file),gsub("_","",nQuire_Ute2$file),nQuire_Ute2$file)
#nQuire_Ute2$file<-test

nQuire_Ute2_flow<-merge(nQuire_Ute2,flow_cytometry,by.x="file",by.y="Sample_ID",all=T)

nQuire_Ute2_flow_sub<-nQuire_Ute2_flow[!(is.na(nQuire_Ute2_flow$Ploidy)|is.na(nQuire_Ute2_flow$Estimated_ploidy)),]

nQuire_Ute2_flow_sub$Estimated_ploidy<-as.numeric(nQuire_Ute2_flow_sub$Estimated_ploidy)

Diff_Ute2<-nQuire_Ute2_flow_sub[nQuire_Ute2_flow_sub$Estimated_ploidy!=nQuire_Ute2_flow_sub$Ploidy,c(1:15,48:49)]
#Kowa001a05
#Kowa001a06
#Zapa002a07

write.xlsx(nQuire_Ute2_flow_sub,"nQuire_flow_cytometry_Ute2.xlsx",overwrite=T)
write.xlsx(nQuire_Ute2,"nQuire_results_Ute2.xlsx",overwrite=T)

