
cat('Welcome!\n
    To use this script;
    \t-First argument must be the mutation location file in bed format. 
      --(You can obtain this by running bedtools intersect -a query.bed -b your.vcf)
    \t-Second argument must be the query bed file.(Usually its peaks are streched).
    \t-Third argument must be the name of this run. This is important for data storage + visualisation
    \t-forth argument must be the path of the output directory. If you wanna run this in the current dir use `pwd`
    ###NOTICE: This script might take some time for highly mutated cancer types. It has not optimised yet.###')
### we need the output of mutation location file. Bed like format done by firats algorithm ###
args=commandArgs()

if (length(args) <8) {
  stop("At least three arguments are required as explained in the intro. Please provide them.", call.=FALSE)
} 


### location of each mutation ###
mutationdosyasi = read.table(file = args[6],sep = '\t',stringsAsFactors = F)

colnames(mutationdosyasi) <- c('chr','start','end','peakno')

#### BED file for querry####
peak_mutation_table=read.table(args[7],sep='\t',header = F,stringsAsFactors = F)
colnames(peak_mutation_table)=c('chr','start','end','peakno')

name=args[8]
output_dir=args[9]

cat('MutationFile=',args[6],'\n')
cat('Query Bed file=',args[7],'\n')
cat('Name of this run=',args[8],'\n')
cat('OutputDir=',args[9],'\n')


mydf=data.frame(peakno='peakno',distance=paste('distance',1:nrow(mutationdosyasi)),
                stringsAsFactors = F) ### this is for getting which difference is belong which peak.
index=0
for ( i in 1:nrow(peak_mutation_table) ){ 
  if (i %% 100 == 0){print(paste(i,':',nrow(peak_mutation_table)))}
  #print(peak_mutation_table[i,"peakno"])
  midpoint <- round( (peak_mutation_table[i,"end"] + peak_mutation_table[i,"start"]+1) / 2)  + 1
  subsetted_mutations <- mutationdosyasi[mutationdosyasi$peakno == peak_mutation_table[i,"peakno"],]
  #print(subsetted_mutations)
  
  if(nrow(subsetted_mutations) != 0) {
    for (j in 1:nrow(subsetted_mutations)){
      index= index+1
      mydf[index,1]=subsetted_mutations[j,'peakno']
      mydf[index,2]=midpoint - subsetted_mutations[j,"start"]
    }
    
  }
}
mydf=as.data.frame(mydf,stringsAsFactors= F)
mydf[,2]=as.numeric(mydf[,2])

### Visualization of the difference between PCA and Other Cancers ###
#PCA=distance_vector
#other = distance_vector
mydf[1]=name

colnames(mydf)=c('typeof','locations')
mydf=mydf[abs(mydf[,"locations"]) < 10000,]

### Density ###
library(ggplot2)
 a <- ggplot(data = mydf,aes(locations,color=typeof,fill=typeof)) +geom_freqpoly(aes(y=..density..),position='identity',binwidth = 75 )+
   xlab(bquote(paste('\nDistance from peak',name,sep = '') ~ (bp))) +ylab(bquote('SNV Mutation Density' ~ (10^-5))) +
   theme(panel.background = element_rect(fill='white'),panel.border = element_blank(),axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
         axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), ### grid options ###
         legend.title = element_blank(),legend.position='top', legend.text = element_text(size = 15),
         axis.text.x=element_text(face = "bold",size = 15) ,
         axis.title.x=element_text(size=20),
         axis.text.y=element_text(face = "bold",size=15),
         axis.title.y=element_text(size=20)) +
   scale_fill_discrete(labels=c('Prostate','Pan Cancer')) + scale_color_discrete(guide=F) +
   scale_y_continuous(labels = function(x)x*100000) + coord_cartesian(xlim=c(-5000,5000))  
   ggsave(plot = a,filename = file.path(output_dir,paste(name,'_snv','_mutation_density','.pdf',sep = '')),height = 105,width = 148,units = 'mm',scale=1.5)

### Histogram ###
# ggplot(data = mydf,aes(locations,color=typeof,fill=typeof)) + geom_histogram(aes(y=..density..),position='identity',alpha=.7,binwidth = 750)  + 
#   xlab(bquote('\nDistance from peak' ~ (bp))) +ylab(bquote('SNV Mutation Density' ~ (10^-5))) +
#   theme(panel.background = element_rect(fill='white'),panel.border = element_blank(),axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
#         axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), ### grid options ###
#         legend.title = element_blank(),legend.position='top', legend.text = element_text(size = 15),
#         axis.text.x=element_text(face = "bold",size = 15) ,
#         axis.title.x=element_text(size=20),
#         axis.text.y=element_text(face = "bold",size=15),
#         axis.title.y=element_text(size=20)) +
#         scale_fill_discrete(labels=c('Prostate','Pan Cancer')) + scale_color_discrete(guide=F) + 
#         scale_y_continuous(labels = function(x)x*100000) + coord_cartesian(xlim=c(-5000,5000))


file_name = paste(name,'locations_dataframe','.txt',sep = '')
write.table(mydf,file = file.path(output_dir,file_name),sep = '\t',quote = F,row.names = F)
save.image(file.path(output_dir,paste(name,'output.RData',sep='_')))
#ggsave('snv_AR_10k_allbinding_TA_Histogram_bin750.pdf',height = 105,width = 148,units = 'mm',scale=1.5)
#save(list = c("PCA","other"),file='AR_10k_tumoronly_mutation.Rdata')


