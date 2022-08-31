#Stampp
options(java.parameters = "-Xmx12g")

require(vcfR)
# convert to genlight 
vcf<-read.vcfR("../Filtered_introgression2/MiasKowaarenosa.vcf", verbose = FALSE)
write.vcf(vcf[vcf@ fix[,1]=="scaffold_1",],file="MiasKowaarenosa_scaffold1.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_2",],file="MiasKowaarenosa_scaffold2.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_3",],file="MiasKowaarenosa_scaffold3.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_4",],file="MiasKowaarenosa_scaffold4.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_5",],file="MiasKowaarenosa_scaffold5.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_6",],file="MiasKowaarenosa_scaffold6.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_7",],file="MiasKowaarenosa_scaffold7.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_8",],file="MiasKowaarenosa_scaffold8.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_9",],file="MiasKowaarenosa_scaffold9.vcf")
