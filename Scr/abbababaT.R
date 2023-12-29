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

#baseDir <- '/home/ntellini/final-ILS-WE'

setwd(baseDir)

# Libraries ----

library(admixr)
library(data.table)

# body ----
Sys.getenv("PATH")


dir.create(path = paste0(baseDir,"/AdmixtoolsVCF/results"),showWarnings = F)

files <- list.files(path = paste0(baseDir,"/AdmixtoolsVCF/"),pattern = ".h.vcf")

for (i in files) {
  
  var <- system(paste0("grep ",i," ",baseDir,"/Sink.txt"),intern = F) 
  
  if ( var == 1 ) {
    
    print(i)
    sink(paste0(baseDir,"/Sink.txt"),append = T)
    print(i)
    sink()
    
    samppref <- paste(strsplit(i,split = "\\.")[[1]][1:2],collapse = ".")
    
    samp <- strsplit(i,split = "\\.")[[1]][1]
    
    pref <- paste0(baseDir,"/AdmixtoolsVCF/",samppref)
    
    snps <- eigenstrat(pref)
    
    d1 <- d(data = snps,W ="WE",X =samp,Y ="SpEU" ,Z = "SJu",params=list(blgsize=0.12))
    
    print(loginfo(d1))
    
    d1 <- as.data.frame(d1)
    
    write.table(x = d1,file = paste0(baseDir,"/AdmixtoolsVCF/results/",samp,".abbababa.res.admixtr"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  }
}
