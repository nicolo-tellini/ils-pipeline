# Tue Jul  6 12:25:21 2021 

# Title: generate chr-by-chr vcfs
# Author: Nicolò T.
# Status: Draft

# Comments:

# Options ----

rm(list = ls())
options(warn = -1)
options(stringsAsFactors = F)
gc()
gcinfo(FALSE)
options(scipen=999)

# Variables ----

ArgsVal <- commandArgs(trailingOnly = T)
baseDir <- ArgsVal[1]

# baseDir <- '/home/nico/pipeline-ILS-RealAssemb'

setwd(baseDir)

# Libraries ----

library(PopGenome)
library(seqinr)
library(data.table)

# body ----

sampdirs <- list.dirs(path = paste0(baseDir,"/AdmixtoolsVCF"),full.names = T,recursive = F)
sampdirs <- sampdirs[grep("results",sampdirs,invert = T)]

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", 
            "chrVII", "chrVIII","chrIX","chrX","chrXI", "chrXII", 
            "chrXIII", "chrXIV", "chrXV", "chrXVI")

for (i in sampdirs) {
  
  samp <- strsplit(i,split = "/")[[1]][length(strsplit(i,split = "/")[[1]])]
  print(samp)
  vcf_file <-list.files(i,pattern = "vcf$")
  vcftab <- data.table::fread(file = paste0(i,"/",vcf_file),sep = "\t")
  
  vcftab <- as.data.frame(vcftab)
  
  colnames(vcftab) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","WE",samp,"SpEU","SJu")
  
  for (j in 1:length(allChr) ) {
  
    tempvcf <- vcftab[vcftab[,1] == j,]
    dir.create(paste0(i,"/",samp,"_",allChr[j]))
    
    vcfr <- list.files(paste0(i,"/",samp,"_",allChr[j]))
    if ( length(vcfr) != 0 ) {
      file.remove(paste0(i,"/",samp,"_",allChr[j],"/",vcfr))
    }
    
    write.table(tempvcf,file = paste0(i,"/",samp,"_",allChr[j],"/",samp,".",allChr[j],".vcf"),append = F,quote = F,sep = "\t",
                col.names = T,row.names = F)
  }
}
