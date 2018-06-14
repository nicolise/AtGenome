#Nicole E Soltis 
#03/28/18
#randomize phenotypes once each to try out permutations for GEMMA

#------------------------------------------------------------------------------
rm(list=ls())

#advice on permutations here:
#https://github.com/genetics-statistics/GEMMA/issues/93

#pipeline note:
#1. run D01_TABtoPEDnMAP.R
#2. copy plink executable to GEMMA_files
#3. in command prompt: cd to GEMMA_files
#4. RUN ./plink --noweb --file D_01_PLINK/dpbinMAF20NA10 --maf 0.2 --make-bed --out binMAF20NA10 } do this ONCE for all AtBcGWAS. NEXT STEP is customized by permutation
#5. run this script (02_FAMaddPhenos.R)
#6. cd to 02_GEMMA/
#7. copy edited _rand.fam, .bim, .bam to 02_GEMMA/
#8. calculate k-matrix with: bash fullrand_GEMMA_kmatrix.sh
#9. run GEMMA: bash fullrand_GEMMA_kmatrix_run.sh

rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/AtGenome/")
#trying with REF lesion phenos for now... where are At transcript phenos?
Phenos <- read.csv("data/00_REF/GWASphenotypesALL.csv")

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("data/01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#first merge Phenos, D.Phenos
#col2 = V2 = accession
names(Phenos)[1] <- "V2"

myFAM_match <- myFAM
Phenos_match <- Phenos

setdiff(myFAM_match$V2,Phenos_match$V2) #none, great!
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#reorder phenos
Phenos_match$V2
myFAM_match$V2

#just keep the first 2 phenotypes
Phenos_match <- Phenos_match[,1:3]

myFAM_match2 <- merge(myFAM_match, Phenos_match, by ="V2")
#also remove null phenotype
myFAM_match2 <- myFAM_match2[,c(1:5,7:ncol(myFAM_match2))]

write.csv(myFAM_match2, "data/02_GEMMA/binMAF20NA10_fam.csv")
write.table(myFAM_match2, "data/02_GEMMA/binMAF20NA10_allphenos.fam", row.names=FALSE, col.names=TRUE)
#copy _allphenos.fam to _rand.fam and also copy 01/.bed and 01/.bim to 02/_rand.
myFAM_match2 <- read.csv("data/02_GEMMA/binMAF20NA10_fam.csv")
