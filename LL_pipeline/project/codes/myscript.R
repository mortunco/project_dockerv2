### aim of this script is to install neccasary dependancies for the files ###
### ggplot2 Installation ###
install.packages("ggplot2",repos="http://cran.rstudio.com/")

library("ggplot2")

source("https://bioconductor.org/biocLite.R")
biocLite("SomaticSignatures")
biocLite("BSgenome.Hsapiens.1000genomes.hs37d5")
