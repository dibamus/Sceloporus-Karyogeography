#This script records phylogenetic signal

library('phytools')
library('vegan')
chr <- read.csv("Data/chromosomes.csv", header = T, stringsAsFactors = F)
tree <- readRDS(file = "Outputs/SceloporusTree.rds")
intersects_df <- readRDS("Outputs/intersects_df.RDS")

tree$tip.label <- paste0(sapply(strsplit(tree$tip.label,"_"), `[`, 1), "_", sapply(strsplit(tree$tip.label,"_"), `[`, 2) )

finalset <- chr[chr$chromosomally_diverse_clade & chr$Chromosome_num != "" & chr$insular_endemic == FALSE & chr$in_GARD,]
finalset$Binomial <- gsub(" ", "_", finalset$Binomial)

tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% finalset$Binomial]);

trait <- finalset$Chromosome_num
names(trait) <- gsub(" ", "_", finalset$Binomial)
trait <- trait[names(trait) %in% tree$tip.label]

# coerce to numeric 
names <- names(trait)
trait <- as.numeric(trait)
names(trait) <- names

trait['Sceloporus uniformis'] <- 26

trait <- trait[match( tree$tip.label, names(trait))]

traitdist <- dist(trait, upper = T, diag = T)
traitdist <- as.matrix(traitdist)

phylodist <- cophenetic.phylo(tree)

# continuous comparison of phylodist v kt dist

mantel(xdis = phylodist, ydis = traitdist, method = "pearson")

plot(phylodist, traitdist, bty = "l", xlab = "phylogenetic distance", ylab = "karyotype distance", pch = 21, bg = "cornflowerblue", cex = 1.5)

## binary comparison of phylodist across kt similarity
group <- c(traitdist)
group[group != 0] <- 1 # if same kt, distance will be 0; if different kt, distance will be >0

par(bty = "l")
boxplot(c(phylodist) ~group , bty = "l", xlab = "karyotype",ylab= "phylogenetic distance", names = c("same", "different"))

aov <- aov(c(phylodist) ~group)

summary(aov)

## phylogenetic signal of kt

 phylosig(tree, trait, method = "K", nsim = 10000, test = T)
 phylosig(tree, trait, method = "lambda", test = T)
