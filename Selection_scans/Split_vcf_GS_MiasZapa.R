#Stampp
options(java.parameters = "-Xmx12g")

require(vcfR)

# convert to genlight 
vcf<-read.vcfR("../Filtered_introgression2/MiasZapaarenosa.vcf", verbose = FALSE)
write.vcf(vcf[vcf@ fix[,1]=="scaffold_1",],file="MiasZapaarenosa_scaffold1.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_2",],file="MiasZapaarenosa_scaffold2.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_3",],file="MiasZapaarenosa_scaffold3.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_4",],file="MiasZapaarenosa_scaffold4.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_5",],file="MiasZapaarenosa_scaffold5.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_6",],file="MiasZapaarenosa_scaffold6.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_7",],file="MiasZapaarenosa_scaffold7.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_8",],file="MiasZapaarenosa_scaffold8.vcf")
write.vcf(vcf[vcf@ fix[,1]=="scaffold_9",],file="MiasZapaarenosa_scaffold9.vcf")

