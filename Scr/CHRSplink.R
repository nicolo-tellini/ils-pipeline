# Tue Jul  6 12:25:21 2021 

# Title:
# Author: Nicolò T.
# Status: Draft

# Comments:

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)
gc()
gcinfo(FALSE)
options(scipen=999)

# Variables ----

ArgsVal <- commandArgs(trailingOnly = T)
baseDir <- ArgsVal[1]

# baseDir <- '/home/ntellini/ILS/pipeline-ILS-WE-WC'

setwd(baseDir)

# Libraries ----

library(data.table)

# body ----

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", 
            "chrVII", "chrVIII","chrIX","chrX","chrXI", "chrXII", 
            "chrXIII", "chrXIV", "chrXV", "chrXVI")

files <- list.files(path = paste0(baseDir,"/AdmixtoolsVCF/"),pattern = "-SpEU-Sj.vcf.RDS")
#done <- read.table(paste0(baseDir,"/todo.txt"))

for (i in files) {
  
  samp <- sapply(strsplit(strsplit(i,split = "\\.")[[1]][1],"-"),"[[",2)
  
# if (sum(done$V1 == samp) == 0) {

  vcf <- readRDS(paste0(baseDir,"/AdmixtoolsVCF/",i))
  
  for ( j in 1:length(allChr) ) {
    print(j)
    vcf[vcf[,1] == allChr[j],1] <- j
  }
  
  fwrite(x = vcf,file = paste0(baseDir,"/AdmixtoolsVCF/",samp,".plink.vcf"),append = F,quote = F,sep = "\t",row.names = F,col.names = T) 
# }
}
