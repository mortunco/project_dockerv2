library(ggplot2)
#### RUN the code below  ###
mydata=read.csv('/project/Batch1/mutation_numbers_per_TF/PeakOverlapWithAR.csv',header = T,stringsAsFactors=F)
random_numbers=read.table(file = '/project/Batch1/final/snv_randomized_mutation_counts.txt',sep = '\t',stringsAsFactors = F)
cat('Random frequency was calculated based on clinical AR raw peaks area and mutations!\n')
clinical_AR_raw_mean=mean(unlist(random_numbers[,4:c(ncol(random_numbers)-1)])/1361367) ## c() around ncol makes this expression working. Other vise it substracts -1 from all of the 4:ncol(random_numbers) numbers. Stupid huh ?
clinical_AR_raw_sd=sd(unlist(random_numbers[,4:c(ncol(random_numbers)-1)])/1361367)
myrandom=rnorm(sd = clinical_AR_raw_sd,mean = clinical_AR_raw_mean,n = 100000)
myrandom=data.frame(values=myrandom,stringsAsFactors = F)
patient_number=102
cat('PatientNumberParameter was set to',patient_number,'Please change it if its not correct.')

#badones <- paste(c('TOP1','GRHL2','MED12','TET2','EP300','RUNX1','H3K9me3',
#                   'H3K4me1','H3K4me3','H3K36me3','H3K27ac'),collapse = '|')

#badones <- paste(c('H3K9me3','H3K4me1','H3K4me3','H3K36me3','H3K27ac','TOP1','motif'),collapse = '|')
#mydata=mydata[!grepl(x = mydata$TF,pattern = badones),]
#library(ggplot2)

### Draws a rectangle that mean + 2sd, mean -2sd to determine outliers. ###
rectangle <- data.frame(xstart=c(clinical_AR_raw_mean + 2*clinical_AR_raw_sd), xend=c(clinical_AR_raw_mean - 2*clinical_AR_raw_sd),col=letters[1]) 
#x$name <- factor(x$name, levels = x$name[order(x$val)])
mydata$TF <- factor(mydata$TF,levels= mydata$TF[order(mydata$MutFreq)])

### Attmept 2 ###
a<-ggplot(data = mydata, aes(y = TF,x=MutFreq,size=MutationFound)) + 
  geom_point(col='dodgerblue2',size=7)  + 
  geom_rect(data = rectangle,inherit.aes = F, aes(xmin = xstart, xmax = xend, 
                                                  ymin = -Inf, ymax = Inf, fill = col), alpha = 0.4,fill='grey34') +
  scale_size_area(breaks=c(0,250,500,1000,2000,3000,6000),max_size = 10)+
  scale_x_continuous(labels = function(x) round(x*1000000/patient_number,2)) +
  #scale_x_continuous(breaks=c(0.1,1.5,2.5))
  labs(x=('Mutation Frequency (per Mb)' ),y='TF')  +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        axis.text.x=element_text(size = 18,color = 'black'),
        axis.title.x=element_text(size=20,face = 'bold',margin=margin(15,0,0,0)),
        axis.text.y=element_text(size=18,color = 'black'),
        axis.title.y=element_text(size=20,face = 'bold',margin=margin(0,15,0,0)) ) 

ggsave(plot=a,file="/project/Batch1/mutation_numbers_per_TF/TF_dotplot_mutation.pdf",height = 8.5,width = 11.7)
save.image(file='/project/Batch1/mutation_numbers_per_TF/myEnvironment.RData')
#write.table(mydata,'/project/Batch1/mutation_numbers_per_TF/tobevisualised.txt',sep='\t',row.names=F,quote=F)
