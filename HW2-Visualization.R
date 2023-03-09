library(tidyverse)

# Open 40 output files
filels <- list.files("output/")

dfls <- gsub(".txt", "", filels) %>% gsub("maize", "m",.) %>%
  gsub("teosinte", "t",.) %>% gsub("increase", "i",.) %>%
  gsub("decrease", "d",.) %>% gsub("-", "",.)

for (i in 1:length(filels)){
  n <- paste("output/", filels[i], sep = "")
  assign(dfls[i], read.table(n, sep = "\t", header = TRUE))
}

# Maize and Teosinte Visualization
## On each chromosome (1-10)
ls.maize <- list(md1, md2, md3, md4, md5, md6, md7, md8, md9, md10)
ls.teosinte <- list(td1, td2, td3, td4, td5, td6, td7, td8, td9, td10)

for (i in 1:10){
  md <- ls.maize[[i]]
  td <- ls.teosinte[[i]]
  m.freq <- md[,4:ncol(md)] %>% unlist() %>% table() %>% as.data.frame()
  m.freq$Group <- "maize"
  t.freq <- td[,4:ncol(td)] %>% unlist() %>% table() %>% as.data.frame()
  t.freq$Group <- "teosinte"
  freq <- rbind(m.freq, t.freq)
  freq$Group <- as.factor(freq$Group)
  n <- paste("Chromosome", i, sep = " ")
  ggplot(freq, aes(x = ., y = Freq, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = n, x = "Genotype", y = "SNP Count")
}

## Across chromosome
maize <- rbind(md1, md2, md3, md4, md5, md6, md7, md8, md9, md10)
teosinte <- rbind(td1, td2, td3, td4, td5, td6, td7, td8, td9, td10)

m.freq <- maize[,4:ncol(maize)] %>% unlist() %>% table() %>% as.data.frame()
m.freq$Group <- "maize"

t.freq <- teosinte[,4:ncol(teosinte)] %>% unlist() %>%
  table() %>% as.data.frame()
t.freq$Group <- "teosinte"

freq <- rbind(m.freq, t.freq)
freq$Group <- as.factor(freq$Group)
ggplot(freq, aes(x = ., y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Across chromosome", x = "Genotype", y = "SNP Count")

# Missing data and amount of heterozygosity
proportion <- list(Group = c(), Genotype = c(), Prop = c())
homozygous <- c("A/A", "C/C", "G/G", "T/T")
missing <- c("-/-")
# Group maize
genotype <- m.freq$.
total <- sum(m.freq$Freq)
# proportion of missing data
miss.idx <- genotype %in% missing
miss.p <- sum(m.freq$Freq[miss.idx]) / total
# proportion of homozygous
homo.idx <- genotype %in% homozygous
homo.p <- sum(m.freq$Freq[homo.idx]) / total
# proportion of heterozygous
hetero.idx <- !(miss.idx | homo.idx)
hetero.p <- sum(m.freq$Freq[hetero.idx]) / total

proportion$Group[1:3] <- rep("maize", 3)
proportion$Genotype[1:3] <- c("missing", "homozygous", "heterozygous")
proportion$Prop[1:3] <- c(miss.p, homo.p, hetero.p)

# Group Teosinte
genotype <- t.freq$.
total <- sum(t.freq$Freq)
# proportion of missing data
miss.idx <- genotype %in% missing
miss.p <- sum(t.freq$Freq[miss.idx]) / total
# proportion of homozygous
homo.idx <- genotype %in% homozygous
homo.p <- sum(t.freq$Freq[homo.idx]) / total
# proportion of heterozygous
hetero.idx <- !(miss.idx | homo.idx)
hetero.p <- sum(t.freq$Freq[hetero.idx]) / total

proportion$Group[4:6] <- rep("teosinte", 3)
proportion$Genotype[4:6] <- c("missing", "homozygous", "heterozygous")
proportion$Prop[4:6] <- c(miss.p, homo.p, hetero.p)

proportion <- as.data.frame(proportion)
proportion$Group <- as.factor(proportion$Group)
proportion$Genotype <- as.factor(proportion$Genotype)
proportion$Prop <- proportion$Prop * 100

ggplot(proportion, aes(x = Genotype, y = Prop, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) + ylim(0, 100) +
  labs(title = "Proportion of Heterozygosity by Group",
       x = "Genotype", y = "%")


