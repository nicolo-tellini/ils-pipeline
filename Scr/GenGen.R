# Mon Jun 21 18:26:17 2021 

# Title: Take SGD genome and create the new aseemblies replacing the alleles with ALT info from variant calling
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

baseDir <- '/home/nico/prog/tempDirSync'

setwd(baseDir)

# Libraries ----

library(seqinr)
library(data.table)

# body ----

dir.create(paste0(baseDir,"/FastaGenomes"),showWarnings = F,recursive = T)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", 
            "chrVII", "chrVIII","chrIX","chrX","chrXI", "chrXII", 
            "chrXIII", "chrXIV", "chrXV", "chrXVI")


allvcfs <- list.files(path = paste0(baseDir,"/VCFs/"),pattern = ".filt.gvcf.recode.vcf.gz$")
allvcfs <- allvcfs[1419:length(allvcfs)]

for (i in allvcfs) {
  
  sgd_genome <- seqinr::read.fasta(file = paste0(baseDir,"/Ref/WE.genome.fa"),set.attributes = F,forceDNAtolower = F)
  
  samp <- sapply(strsplit(i,split = "\\."),"[[",1)
  
  vcf <- fread(file = paste0(baseDir,"/VCFs/",i))
  vcf <- as.data.frame(vcf)
  
  vcf <- vcf[as.numeric(vcf[,6]) > 20,]
  
  vcf <- vcf[grep("0/0",vcf[,10],invert = T),]
  
  vcf$counts <- apply(as.data.frame(vcf[,5]),1,function(x) length(strsplit(x,"")[[1]]))
  vcf <- vcf[vcf$counts == 1,]
  vcf$counts <- apply(as.data.frame(vcf[,4]),1,function(x) length(strsplit(x,"")[[1]]))
  vcf <- vcf[vcf$counts == 1,]
  
  vcf[,10] <- sapply(strsplit(vcf[,10] ,":"),"[[",1)
  
  vcf <- vcf[grep(",",vcf[,4],invert = T),]
  vcf <- vcf[grep(",",vcf[,5],invert = T),]
  
  for (j in 1:nrow(vcf)) {
    sgd_genome[[vcf[j,1]]][vcf[j,2]] <- vcf[j,5]
  }
  
  seqinr::write.fasta(sequences = sgd_genome,names = names(sgd_genome),open = "w",file.out = paste0(baseDir,"/FastaGenomes/",samp,".genome.fa"))
}

