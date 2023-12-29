# Wed Jun 30 19:23:47 2021 

# Title: VCF for ADMIXTOOLS
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
i <- ArgsVal[2]

# baseDir <- '/home/ntellini/ILS/pipeline-ILS-WE-WC'

setwd(baseDir)

# Libraries ----

library(seqinr)
library(data.table)
library(stringr)

# Body ----

#dir.create(paste0(baseDir,"/AdmixtoolsVCF"),showWarnings = F,recursive = T)

#allsnps <- list.files(path = paste0(baseDir,"/SNP"),pattern = ".snps$",full.names = T)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", 
            "chrVII", "chrVIII","chrIX","chrX","chrXI", "chrXII", 
            "chrXIII", "chrXIV", "chrXV", "chrXVI")

#done <- read.table(paste0(baseDir,"/already_done.txt"))

#for (i in allsnps) {
  
  print(i) 
  lenp <- length(strsplit(i,split = "/")[[1]])
  samp <- strsplit(strsplit(i,split = "/")[[1]][lenp],split = "-")[[1]][2]
  
#if (sum(done$V1 == samp) == 0) {
    
    snp <- fread(file = i,sep = "\t",header = T,fill = T,data.table = FALSE,nThread = 8)
    
    snp <- snp[grep("-",snp[,1],invert = T),]
    
    if ( nrow(snp) > 0) {
      
      snp$count <- nchar(snp[,1])
      snp <- snp[snp[,"count"] == 4,]
      snp$count <- NULL
      
      if ( nrow(snp) > 0){
        
        snp$pos1 <- substr(snp[,1], 1, 1)
        snp$pos2 <- substr(snp[,1], 2, 2)
        snp$pos3 <- substr(snp[,1], 3, 3)
        snp$pos4 <- substr(snp[,1], 4, 4)
        snp$pos4 <- toupper(snp$pos4)
        
        # keep only biallelic ----
        
        snp$counts <- 0
        snp$counts <- apply(snp[,14:17],1,function(x) length(unique(x)) )
        print(nrow(snp[snp$counts != 2,]))
        snp <- snp[snp$counts == 2,]
        
        # keep only those on the same chromosome exept for S. jurei because it is rearranged ----
        if ( nrow(snp) > 0){
          
          snp <- snp[snp$sequence_1_Contig == snp$sequence_2_Contig & snp$sequence_2_Contig == snp$sequence_3_Contig,]
          
          # remove variants on WE telomere ----
          
          cbsubtel <- read.table(file = paste0(baseDir,"/Rep/Ann/CBS432.subtel.txt"))
          
          cbsubtel$chrs <- allChr
          
          vcfsnp <- data.frame()
          for (k in 1:nrow(cbsubtel)) {
            print(k)
            Tsnp <- snp[snp$sequence_3_Contig == cbsubtel[k,"chrs"] & snp$sequence_3_PosInContg >= cbsubtel[k,1] & snp$sequence_3_PosInContg <= cbsubtel[k,2],]
            vcfsnp <- rbind(vcfsnp,Tsnp)
          }
          rm(snp,Tsnp)
          
          if (nrow(vcfsnp) > 0){
            
            # vcfsnp$alleles <- as.character(apply(vcfsnp[,14:17],1,function(x) paste(unique(x),collapse = ",")))
            # 
            # vcfsnp <-  vcfsnp[!is.na(vcfsnp$alleles),]
            # 
            # vcfsnp <-  vcfsnp[vcfsnp$alleles != "NA",]
            
            # CHR P2
            CHR <- vcfsnp[,5]
            # POS P2
            POS <- vcfsnp[,6]
            # ID P2
            ID <- rep(".",times=nrow(vcfsnp))
            # REF P2
            REF <- vcfsnp$pos2
            
            P1booleP2 <- vcfsnp$pos1 == vcfsnp$pos2
            P2booleP2 <- vcfsnp$pos2 == vcfsnp$pos2
            P3booleP2 <- vcfsnp$pos3 == vcfsnp$pos2
            P4booleP2 <- vcfsnp$pos4 == vcfsnp$pos2
            
            P1 <- as.numeric(!P1booleP2)
            P2 <- as.numeric(!P2booleP2)
            P3 <- as.numeric(!P3booleP2)
            P4 <- as.numeric(!P4booleP2)
            
            vcfsnp[vcfsnp$pos1 == vcfsnp$pos2,"pos1"] <- "0"
            vcfsnp[vcfsnp$pos3 == vcfsnp$pos2,"pos3"] <- "0"
            vcfsnp[vcfsnp$pos4 == vcfsnp$pos2,"pos4"] <- "0"
            vcfsnp[,"pos2"] <- "0"
            
            ALT <- rep("",times=nrow(vcfsnp))
            
            for (l in 14:17) {
              
              rows <- which(!grepl("0", as.matrix(vcfsnp[,l]), fixed=TRUE))
              ALT[rows] <- vcfsnp[rows,l]
              
            }
          
            QUAL <- rep(100,times=nrow(vcfsnp)) 
            FILTER <- rep(".",times=nrow(vcfsnp))  
            INFO  <- rep(".",times=nrow(vcfsnp))
            FORMAT <- rep("GT",times=nrow(vcfsnp))
            
            vcf <- data.frame(CHR,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,P1,P2,P3,P4)
            
            # rename columns ----
            
            names(vcf)[1] <- "#CHROM"
            colnames(vcf)[10:13] <- c("WE",samp,"SpEU","SJu")
            
            saveRDS(vcf, file = paste0(baseDir,"/AdmixtoolsVCF/","WE-",samp,"-SpEU-Sj.vcf.RDS")) 
            
          }
        }
      }
    }
#  }
#}
