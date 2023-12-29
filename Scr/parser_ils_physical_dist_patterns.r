# Fri Oct  1 17:38:03 2021 

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

baseDir <- '/home/nico/prog/tempDirSync/SNP'
# allDir <- list.dirs()
# allFiles <- list.files()
setwd(baseDir)

# Libraries ----

library(data.table)
library(ggplot2)
library(Biostrings)
library(GenomicRanges)
library(ggpubr)
library(plotly)

# body ----

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", 
            "chrVII", "chrVIII","chrIX","chrX","chrXI", "chrXII", 
            "chrXIII", "chrXIV", "chrXV", "chrXVI")

cbsubtel <- read.table(file = "/home/nico/prog/tempDirSync/Rep/Ann/CBS432.subtel.txt")

cbsubtel$chrs <- allChr

all_files <- list.files()

dir.create(path = "/home/nico/prog/tempDirSync/WG_patterns")

for (f in all_files) {
  
  snp <- fread(file = f,sep = "\t",header = T,fill = T)
  snp <- as.data.frame(snp)
  
  snp <- snp[grep("-",snp[,1],invert = T),]
  
  snp$count <- lapply(strsplit(snp[,1],split = ""),FUN = function(x)length(x))
  snp <- snp[snp[,"count"] == 4,]
  snp$count <- NULL
  
  snp$pos1 <- sapply(strsplit(snp[,1],split = ""),"[[",1)
  snp$pos2 <- sapply(strsplit(snp[,1],split = ""),"[[",2)
  snp$pos3 <- sapply(strsplit(snp[,1],split = ""),"[[",3)
  snp$pos4 <- toupper(sapply(strsplit(snp[,1],split = ""),"[[",4))
  
  # keep only biallelic ----
  
  snp$counts <- apply(snp[,14:17],1,function(x) length(unique(x)))
  print(nrow(snp[snp$counts != 2,]))
  snp <- snp[snp$counts == 2,]
  
  # keep only those on the same chromosome exept for S. jurei because it is rearranged ----
  
  snp <- snp[snp$sequence_1_Contig == snp$sequence_2_Contig & snp$sequence_2_Contig == snp$sequence_3_Contig,]
  
  # remove variants on WE telomere ----
  
  vcfsnp <- data.frame()
  for (k in 1:nrow(cbsubtel)) {
    print(k)
    Tsnp <- snp[snp$sequence_3_Contig == cbsubtel[k,"chrs"] & snp$sequence_3_PosInContg >= cbsubtel[k,1] &  snp$sequence_3_PosInContg <= cbsubtel[k,2],]
    vcfsnp <- rbind(vcfsnp,Tsnp)
  }
  rm(snp,Tsnp)
  
  A <- vcfsnp$pos4
  
  vcfsnp[vcfsnp$pos1 == A,"pos1"] <- "A"
  vcfsnp[vcfsnp$pos1 != A & vcfsnp$pos1 != "A" ,"pos1"] <- "B"
  
  vcfsnp[vcfsnp$pos2 == A,"pos2"] <- "A"
  vcfsnp[vcfsnp$pos2 != A & vcfsnp$pos2 != "A" ,"pos2"] <- "B"
  
  vcfsnp[vcfsnp$pos3 == A,"pos3"] <- "A"
  vcfsnp[vcfsnp$pos3 != A & vcfsnp$pos3 != "A" ,"pos3"] <- "B"
  
  vcfsnp$pos4 <- "A"
  
  vcfsnp$pattern <- apply(vcfsnp[,14:17],1,function(x)paste(x,sep = "",collapse = ""))  
  
  BBBA <- vcfsnp[vcfsnp$pattern == "BBBA",]
  AABA <- vcfsnp[vcfsnp$pattern == "AABA",]  
  AAAA <- vcfsnp[vcfsnp$pattern == "AAAA",]  
  ABBA <- vcfsnp[vcfsnp$pattern == "ABBA",]  
  BABA <- vcfsnp[vcfsnp$pattern == "BABA",]
  ABAA <- vcfsnp[vcfsnp$pattern == "ABAA",]  
  BAAA <- vcfsnp[vcfsnp$pattern == "BAAA",]  

  mygenome <- readDNAStringSet("/home/nico/prog/tempDirSync/Ref/Asm/CBS432.genome.fa")
  chrSizes <- width(mygenome)
  names(chrSizes) <- names(mygenome)
  print(chrSizes)
  
  bins   <- as.data.frame(tileGenome(chrSizes, tilewidth=1000, cut.last.tile.in.chrom=T))
  
  bins$BABA <- ""
  bins$BBBA <- ""
  bins$AABA <- ""
  bins$AAAA <- ""
  bins$ABBA <- ""
  bins$ABAA <- ""
  bins$BAAA <- ""
  
  for (i in 1:nrow(bins)) {
    print(i)
    
    bins[i,"BBBA"] <- length(BBBA[BBBA[,8]==bins[i,1] & BBBA[,9]>= bins[i,2] & BBBA[,9] <= bins[i,3],1])
    bins[i,"AABA"] <- length(AABA[AABA[,8]==bins[i,1] & AABA[,9]>= bins[i,2] & AABA[,9] <= bins[i,3],1])
    bins[i,"AAAA"] <- length(AAAA[AAAA[,8]==bins[i,1] & AAAA[,9]>= bins[i,2] & AAAA[,9] <= bins[i,3],1])
    bins[i,"BABA"] <- length(BABA[BABA[,8]==bins[i,1] & BABA[,9]>= bins[i,2] & BABA[,9] <= bins[i,3],1])
    bins[i,"ABBA"] <- length(ABBA[ABBA[,8]==bins[i,1] & ABBA[,9]>= bins[i,2] & ABBA[,9] <= bins[i,3],1])
    bins[i,"ABAA"] <- length(ABAA[ABAA[,8]==bins[i,1] & ABAA[,9]>= bins[i,2] & ABAA[,9] <= bins[i,3],1])
    bins[i,"BAAA"] <- length(BAAA[BAAA[,8]==bins[i,1] & BAAA[,9]>= bins[i,2] & BAAA[,9] <= bins[i,3],1])
    
  }
  
  strain <- base::strsplit(x = f,split = "-")[[1]][2]
  
  write.table(bins,file = paste0("/home/nico/prog/tempDirSync/WG_patterns/",strain,".CBS432coord.bins.patterns"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  
  for (j in 6:12) {
    bins[,j] <- as.numeric(bins[,j])
    bins[,j] <- round(as.numeric(bins[,j]) / max(bins[,j]),digits = 2)
  }
  
  bins$strand <- (as.numeric(bins$end) + as.numeric(bins$start)) / 2 
  
  for (p in allChr) {
    
    bins2 <- bins[bins$seqnames == p,]
    
    a <- ggplot() + 
      geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"BABA"])) +
      ylim(c(0,1)) + xlab("") + ylab("num. BABA sites /\nmax val BABA") +
      theme_bw()     
    b <-  ggplot() + 
      geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"ABBA"])) +
      ylim(c(0,1)) + xlab("") + ylab("num. ABBA sites /\nmax val ABBA") +
      theme_bw()
    c <-   ggplot() + 
      geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"BAAA"])) +
      ylim(c(0,1)) + xlab("") + ylab("num. BAAA sites /\nmax val BAAA") +
      theme_bw()
    d <-  ggplot() + 
      geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"ABAA"])) +
      ylim(c(0,1)) + xlab("") + ylab("num. ABAA sites /\nmax val ABAA") +
      theme_bw()
    e <-  ggplot() + 
      geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"AABA"])) +
      ylim(c(0,1)) + xlab("") + ylab("num. AABA sites /\nmax val AABA") +
      theme_bw()
    
    f <-  ggplot() + 
      geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"AAAA"])) +
      ylim(c(0,1)) + xlab("") + ylab("num. AAAA sites /\nmax val AAAA") +
      theme_bw()
    
    g <-  ggplot() + 
      geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"BBBA"])) +
      ylim(c(0,1)) + xlab(paste0("1kb non-over. windonws on ",p, " of ", strain) ) + ylab("num. BBBA sites /\nmax val BBBA") +
      theme_bw() 
    
    # p2 <- ggplot() +
    #   geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"BABA"]),color="blue") +
    #   geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"ABBA"]),color="black") +
    #   geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"BAAA"]),color="orange") +
    #   geom_line(bins2,mapping = aes(x=bins2[,"strand"],y=bins2[,"ABAA"]),color="red") +
    #   ylim(c(0,1)) + xlab(paste0("1kb non-over. windonws on ",p, " of ", strain) ) +
    #   theme_bw()
    
    plot <- ggarrange(a,b,c,d,e,f,g,nrow = 7)
    # html <- ggplotly(a)
    #html <- ggplotly(p2)
    
    pPath <- file.path(paste0("/home/nico/prog/tempDirSync/WG_patterns/",strain,".",p,".sites.pdf"))
    pdf(file = pPath, width = 15, height = 10)
    print(plot)
    dev.off()
    
    # htmlp <- ggplotly(p2)
    # htmlwidgets::saveWidget(as_widget(htmlp),paste0("/home/nico/prog/tempDirSync/WG_patterns/",strain,".",p,".sites.pdf"))
  }
}

allChr <- c()