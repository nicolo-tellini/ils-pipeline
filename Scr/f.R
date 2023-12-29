  # Tue Jul  6 12:25:21 2021 
  
  # Title:
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
  
 # baseDir <- '/home/nico/pipeline-ILS'
  
  setwd(baseDir)
  
  # Libraries ----
  
  library(PopGenome)
  library(seqinr)
  
  # body ----
  
  sampdirs <- list.dirs(path = paste0(baseDir,"/AdmixtoolsVCF"),full.names = T,recursive = F)
  sampdirs <- sampdirs[grep("results",sampdirs,invert = T)]
  
  allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", 
              "chrVII", "chrVIII","chrIX","chrX","chrXI", "chrXII", 
              "chrXIII", "chrXIV", "chrXV", "chrXVI")
  
  sampdirs <- paste0(sampdirs)
  
  for (i in sampdirs) {
    
    samp <- strsplit(i,split = "/")[[1]][length(strsplit(i,split = "/")[[1]])]
  
    txtcheck <- list.files(i,pattern = "txt")
    
    if ( length(txtcheck) != 0) {
      file.remove(paste0(i,"/",txtcheck))
    }
    
    txtcheck <- list.files(i,pattern = "txt")
    
    if ( length(txtcheck) == 0 ) {
      
      print(samp)
      
      vcf <- readData(path = i, format = "VCF",SNP.DATA = T,big.data=T)
      
      vcf <- set.populations(vcf,list("WE",samp,"SpEU"),diploid = T)
      
      vcf <- set.outgroup(vcf,"SJu",diploid = T)
      
      ### do Martin's f statistic ----
      
      vcf <- introgression.stats(vcf,do.D = T, do.df=F,do.RNDmin = F, block.size = F, keep.site.info = TRUE)
      
      # D value
      pig <- vcf@f * 100 # percentage introgressed genome 
      
      df <- data.frame(samp,pig)
      
      colnames(df) <- c("sample","f")
      
      write.table(df,file = paste0(i,"/",samp,".f.txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)
      
      print(paste0(samp," f genereated"))
    }
  }