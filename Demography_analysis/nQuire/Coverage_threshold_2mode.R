
filelist = list.files(pattern = "*.statistics$") # Make a file list from all count text files

for (file in filelist)
	{assign(file,as.data.frame(t(read.table(file,sep=",",header=F))))
	assign(file,get(file)[-1,])
	}

index_max<-c()

for (file in filelist)
	{file_now<-get(file)
	file_now$V2<-as.numeric(file_now$V2)
	diff=c()
	for (i in 2:length(file_now$V2))
	{diff<-c(diff,file_now$V2[i]-file_now$V2[i-1])
	}
	max1<-min(which(diff[min(which(diff>0)):length(diff)]<0))+min(which(diff>0))
	index_max<-c(index_max,max1*2)
	}

write.table(data.frame(cbind(filelist,index_max)),"Coverage_threshold_2mode.table",row.names=F,sep="\t",quote=F)
