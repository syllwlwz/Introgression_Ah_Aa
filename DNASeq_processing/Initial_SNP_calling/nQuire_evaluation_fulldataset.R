require(openxlsx)
Species<-read.xlsx("Sequencing_list_introgression_formerge.xlsx",1)
nQuire<-read.table("nQuire_depth_summary_diff.txt")
colnames(nQuire)<-c("Sample_ID","Percentage_SNPs_used","Min_coverage_used","Max_coverage_used","nQuire_ploidy","Flow_cytometry_ploidy","Unmatching_ploidy")
all<-merge(nQuire,Species,by="Sample_ID",all.x=T)
write.xlsx(all,"nQuire_evaluation.xlsx")






