setwd('~/project/Batch1/mutation_numbers_per_TF/')
myfiles <- list.files(pattern = '.txt')
#print(myfiles)
#print(myfiles)
mydf = data.frame(AreaCovered=numeric(),MedianPeak=numeric(),MutationNumber=numeric())
mynames=c()
for (i in 1:length(myfiles)){
print(myfiles[i])
temp=read.table(myfiles[i],sep = '\t',stringsAsFactors = F)
AreaCovered= sum(temp$V3 - temp$V2)
MedianPeak=median(temp$V3 - temp$V2)
MutationNumber=sum(temp[,ncol(temp)])
tfname <- myfiles[i]
mynames=c(mynames,tfname)
mytempdf <- data.frame(AreaCovered=AreaCovered,MedianPeak=MedianPeak,MutationNumber=MutationNumber,stringsAsFactors=F)
mydf=rbind.data.frame(mydf,mytempdf,stringsAsFactors=F)
}
colnames(mydf) <- c('AreaCovered','MedianPeak',"MutationNumber")
mydf$MutFreq <- mydf$MutationNumber / mydf$AreaCovered
mydf$TF <- mynames
write.table(x = mydf,file = 'PeakOverlapWithAR.csv',row.names = F,quote = F,col.names=T,sep=',')
