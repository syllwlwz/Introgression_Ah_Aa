
filelist = list.files(pattern = "*.statistics$") # Make a file list from all count text files

for (file in filelist)
	{assign(file,as.data.frame(t(read.table(file,sep=",",header=F))))
	assign(file,get(file)[-1,])
	}

index_max<-c()
for (file in filelist)
	{file_now<-get(file)
	file_now$V2<-as.numeric(file_now$V2)
	index_max<-c(index_max,min(which(file_now$V2<100000))-1)
	}

write.table(data.frame(cbind(filelist,index_max)),"Coverage_threshold.table",row.names=F,sep="\t",quote=F)
