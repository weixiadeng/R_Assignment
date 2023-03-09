library(tidyverse)
# Read fang_et_al_genotypes file
fang <- read.table("fang_et_al_genotypes.txt", sep = "\t", header = TRUE)
# Group of maize and teosinte
maize <- c("ZMMIL", "ZMMLR", "ZMMMR")
teosinte <- c("ZMPBA", "ZMPIL", "ZMPJA")
# Find index of maize and teosinte from third column of fang_et_al_genotypes
maize.idx <- fang$Group %in% maize
teosinte.idx <- fang$Group %in% teosinte
# Subset maize and teosinte from fang_et_al_genotypes data
# Only keep genotypes data for each group 
# (drop Sample_ID, JG_OTU and Group columns)
fang.maize <- fang[maize.idx, -c(1:3)]
fang.teosinte <- fang[teosinte.idx, -c(1:3)]
# Add header of maize and teosinte (genotypes) to Subset of data
fang.maize <- rbind(colnames(fang.maize), fang.maize)
fang.teosinte <- rbind(colnames(fang.teosinte), fang.teosinte)
# Transpose two subsets of genotype data
# First column is the genotype
maize.t <- t(fang.maize) %>% as.data.frame()
teosinte.t <- t(fang.teosinte) %>% as.data.frame()

# Read snp_position file
snp <- read.table("snp_position.txt", sep = "\t", header = TRUE)
# Only keep SNP id (first column),
# chromosome location (third column),
# nucleotide location (fourth column)
snp.sub <- snp[,c(1,3,4)]
# Remove Position values are "", "multiple", "unknown"
position.idx <- snp.sub$Position %in% c("", "multiple", "unknown")
snp.sub <- snp.sub[!position.idx,]
# Set chromosome as factor, Position as numeric
snp.sub$Chromosome <- as.factor(snp.sub$Chromosome)
snp.sub$Position <- as.numeric(snp.sub$Position)

# Maize data
for (i in 1:10){
  # Subset SNP by Chromosome 1 to 10
  chromosome.idx <- snp.sub$Chromosome == i
  # Merge subset SNP data with maize genotype data by genotype
  df <- merge(x = snp.sub[chromosome.idx,], y = maize.t,
              by.x = "SNP_ID", by.y = "1")
  # Sort position by increasing order
  df.1 <- df[order(df$Position, decreasing = FALSE),]
  # Save df.1 to output folder
  n.1 <- paste("output/maize-increase-", i, ".txt", sep = "")
  write.table(df.1, file = n.1, sep = "\t")
  # Sort position by decreasing order
  df.2 <- df[order(df$Position, decreasing = TRUE),]
  # Replace missing data "?" to "-"
  df.2[,4:ncol(df.2)] <- lapply(df.2[,4:ncol(df.2)],
                                function(x) str_replace_all(x, "\\?", "-"))
  # Save df.2 to output folder
  n.2 <- paste("output/maize-decrease-", i, ".txt", sep = "")
  write.table(df.2, file = n.2, sep = "\t")
}

# Teosinte data
for (i in 1:10){
  # Subset SNP by Chromosome 1 to 10
  chromosome.idx <- snp.sub$Chromosome == i
  # Merge subset SNP data with teosinte genotype data by genotype
  df <- merge(x = snp.sub[chromosome.idx,], y = teosinte.t,
              by.x = "SNP_ID", by.y = "1")
  # Sort position by increasing order
  df.1 <- df[order(df$Position, decreasing = FALSE),]
  # Save df.1 to output folder
  n.1 <- paste("output/teosinte-increase-", i, ".txt", sep = "")
  write.table(df.1, file = n.1, sep = "\t")
  # Sort position by decreasing order
  df.2 <- df[order(df$Position, decreasing = TRUE),]
  # Replace missing data "?" to "-"
  df.2[,4:ncol(df.2)] <- lapply(df.2[,4:ncol(df.2)],
                                function(x) str_replace_all(x, "\\?", "-"))
  # Save df.2 to output folder
  n.2 <- paste("output/teosinte-decrease-", i, ".txt", sep = "")
  write.table(df.2, file = n.2, sep = "\t")
}
