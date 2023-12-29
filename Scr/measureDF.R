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

#baseDir <- '/home/nico/pipeline-ILS-RealAssemb'

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

chrlen <- read.table(paste0(baseDir,"/Rep/Ann/WE.chrs.txt"))

for (i in sampdirs) {
  
  print(i)
  samp <- strsplit(i,split = "/")[[1]][length(strsplit(i,split = "/")[[1]])]
  
  dirs <- list.dirs(i,full.names = T)[grep("chr",list.dirs(i,full.names = T))]
  
  for (p in dirs) {
    
    pdfcheck <- list.files(p,pattern = "pdf")
    
    if ( length(pdfcheck) != 0) {
      file.remove(paste0(p,"/",pdfcheck))
    }
    
    pdfcheck <- list.files(p,pattern = "pdf")
    
    if ( length(pdfcheck) == 0 ) {
      
      chr <- strsplit(strsplit(p,split = "/")[[1]][length(strsplit(p,split = "/")[[1]])],split = "_")[[1]][2]
      print(chr)
      vcf <- readData(path = p, format = "VCF",SNP.DATA = T)
      
      # if ( vcf@n.biallelic.sites > 0 ) {
      
      vcf <- set.populations(vcf,list("WE",samp,"SpEU"),diploid = T)
      
      vcf <- set.outgroup(vcf,"SJu",diploid = T)
      
      ### do.df		Bd-fraction ----
      
      # vcf$df	[3] Bastian Pfeifer and Durrell D. Kapan (2019). 
      # Estimates of introgression as a function of pairwise distances. 
      # BMC Bioinformatics. https://doi.org/10.1186/s12859-019-2747-z
      
      vcf <- sliding.window.transform(vcf,jump = 1000,width=1000,start.pos=1,end.pos = chrlen[chrlen[,1] == chr,2],type = 2)
      
      vcf <- introgression.stats(vcf, l.smooth=T,do.df = T)
      
      vcf <- weighted.jackknife(vcf,per.region = F)
      
      ind <- which( abs(vcf@df) >= 5 * mean( abs(vcf@df[!is.na(vcf@df)]) ) ) # red dots 
      
      # ind2 <- which((1-adjP) < 0.01) 
      
      pdfname <- paste0(p,"/",samp,".",chr,".dfplot.pdf")
      
      write.table(x = vcf@df,append = F,file = paste0(p,"/",samp,".",chr,".df_values.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

      # Open a pdf file
      pdf(pdfname) 
      # 2. Create a plot
      plot(1:length(vcf@df),as.numeric(vcf@df),type="l",ylim=c(-1,1),xlab=c("1 Kb non overlapping windows"),ylab=c("df"),
           col="blue",lwd=2,main= paste0(samp,".",chr,".dfplot.pdf"))
      if (length(ind) > 1) {
        len <- 1:length(vcf@df)
        y <- (as.numeric(vcf@df[ind]))
        points(len[ind],y,col="red",cex=0.5,pch=16)
      }
      # Close the pdf file
      dev.off() 
    }
  }
}
