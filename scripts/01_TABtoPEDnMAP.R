#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/AtGenome/")
library(tidyr)
#convert each file to binary
SNPsMAF20 <- read.csv("data/JulinsAcc3_MAF20_RFF.csv")
mySNPs <- SNPsMAF20

#remove all SNPs with NA in > 9 accessions (10%)
mySNPs$NAcount <- apply(mySNPs[,c(2:96)], 1, function(x) sum(is.na(x)))
mySNPs <- mySNPs[mySNPs$NAcount <= 9,] #drops none
unique(mySNPs$NAcount) #none!

#and now for making PED format for PLINK!
  #do not need positional info: just SNP states for PED
#turn df sideways (individuals as rows, SNPs as columns)
#split each genotype into 2 identical columns (PED assumes diploid)
#add a first column: FAM1 (no info on isolate families)
#second column: isolate ID
#third column: father ID (a column of zeros)
#fourth column: mother ID (a column of zeros)
#fifth column: individual sex = 1 (all assumed same)
#sixth  column: binary  phenotype (all = 1)
#fix column order
mySNPs2 <- mySNPs[,-c(1,97)]

#recode 0 as 2 so not "missing" SNP
mySNPs2[mySNPs2==0]<- 2

#turn all SNPs to "diploid"
#haha, it takes 4 days to do this as a "for" loop (for each row, rbind twice)
#because is.na <-0 before this step, there should be NO heterozygous SNP calls
#this is super fast:
mySNPs3 <- mySNPs2[rep(1:nrow(mySNPs2),each=2),] 

#transpose and format for PED
mySNPs4 <- as.data.frame(t(mySNPs3))
#add binary phenotype = 1 (6)
mySNPs4 <- cbind("Pheno" = 1, mySNPs4)
#add individual sex = 1 (5)
mySNPs4 <- cbind("sex" = 1, mySNPs4)
#add Mother = 0 (4)
mySNPs4 <- cbind("Mother" = 0, mySNPs4)
#add Father = 0 (3)
mySNPs4 <- cbind("Father" = 0, mySNPs4)
#turn row names into column 2
mySNPs4 <- cbind(rownames(mySNPs4), mySNPs4)
colnames(mySNPs4)[1] <- 'Genotype'
#add the fam column (1)
mySNPs4 <- cbind("FAM" = "FAM1", mySNPs4)
myPED <- mySNPs4


#add a phenotype for PED? 
#NA is fine for missing phenotypes
#since many phenotypes, just add as consecutive columns to *.fam, and run GEMMA in a loop over phenotypes

#make a MAP file for plink (need it to make the bed (binary ped) file from ped)
myMAP <- as.data.frame(mySNPs[,c("SNP")])
#split SNP into chrom and pos
names(myMAP)[1] <- "SNP"
myMAP2 <- separate(data = myMAP, col = SNP, into = c("Chr", "Pos"))
myMAP2$SNPID <- paste("SNP",myMAP2$Pos, sep="")
myMAP2$SNPcM <- 0
myMAP2 <- myMAP2[,c(1,3,4,2)]
write.table(myMAP2, "data/01_PLINK/dpbinMAF20NA10.map", row.names=FALSE, col.names=FALSE)
write.table(myPED, "data/01_PLINK/dpbinMAF20NA10.ped", row.names=FALSE, col.names=FALSE)

