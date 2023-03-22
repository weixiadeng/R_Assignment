library(tidyverse)

# genotype data
file.info("fang_et_al_genotypes.txt")$size
fang <- read.table("fang_et_al_genotypes.txt", sep = "\t", header = TRUE)
dim(fang)
fang$Group <- as.factor(fang$Group)
table(fang$Group)
str(fang[,1:15])

# snp data
file.info("snp_position.txt")$size
snp <- read.table("snp_position.txt", sep = "\t", header = TRUE)
dim(snp)
str(snp)
