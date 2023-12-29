# Fri Jan 14 10:39:31 2022 

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

baseDir <- '/home/nico/temp/AdmixtoolsVCF'
# allDir <- list.dirs()
# allFiles <- list.files()
setwd(baseDir)

# Libraries ----

library(dplyr)
library(readr)

# body ----

vcflist <- list.files(pattern = "-SpEU-Sj.vcf")

for (i in vcflist) {
  
  vcf <- data.table::fread(sep = "\t", cmd = paste("cut -f1,2,4,5,10,11,12,13 ",i))
  vcf <- as.data.frame(vcf)
  vcf$pattern <- apply(vcf[,5:8],1,function(x)(paste(x,sep = "",collapse = "")))
  
  str <- strsplit(i,split = "-")[[1]][2]
  
  target <- c("0110","1001")
  vcf %>%
    filter(pattern %in% target ) %>%
    write_tsv(paste0("abba-sites-",str,".tsv"))
  
  target <- c("0101","1010")
  vcf %>%
    filter(pattern %in% target) %>%
    write_tsv(paste0("baba-sites-",str,".tsv")) 
}
