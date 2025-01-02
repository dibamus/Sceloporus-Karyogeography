#### SPECIES-LEVEL ANALYSIS ####
library("sf")
library("tidyverse")
SceloporusRanges <- readRDS("Outputs/Sceloporus_data_frame.RDS")
# Which species' ranges intersect?
# Which Species Ranges intersect (According to GARD data)
ScelIntersects <- st_intersects(SceloporusRanges, sparse = FALSE)
row.names(ScelIntersects)<- gsub(" ","_",SceloporusRanges$species)
colnames(ScelIntersects)<- gsub(" ","_",SceloporusRanges$species)

#Which species' ranges intersect according to Rivera et al?
Riv <- read.csv("Data/Rivera species overlaps.csv")

Riv_matrix <- sapply(1:length(Riv$Focalspecies), function(x){
  sapply(Riv$Focalspecies, function(y){
    y %in% Riv[x,2:45]
  })
} )

colnames(Riv_matrix) <- row.names(Riv_matrix) 

# Which species have the same karyotypes?
data <- SceloporusRanges$chromosome_group

f <- function(x) x == data

KarIntersects <- sapply(data, f, simplify = "array")
colnames(KarIntersects) <- SceloporusRanges$species
rownames(KarIntersects) <- SceloporusRanges$species

#Make a big df of all species pairs and where they intersect

mx <- lower.tri(ScelIntersects)
row.names(mx) <- row.names(ScelIntersects)
colnames(mx) <- colnames(ScelIntersects)

intersects_df <- as.data.frame(which(mx, arr.ind = TRUE, useNames = TRUE))
intersects_df$row <- colnames(ScelIntersects)[intersects_df$row]
intersects_df$col <- colnames(ScelIntersects)[intersects_df$col]
intersects_df$Karyotype <- KarIntersects[lower.tri(KarIntersects, diag = FALSE)]
intersects_df$GARD <- ScelIntersects[lower.tri(ScelIntersects, diag = FALSE)]

intersects_df$Rivera <- sapply(1:dim(intersects_df)[1], function(x){
  if(all(c(intersects_df$row[x],intersects_df$col[x]) %in% colnames(Riv_matrix))){
    return(Riv_matrix[intersects_df$row[x],intersects_df$col[x]])
  }
  else{return(NA)}
  
})
#for any list of species names, mark the rows in which both those species are included
relevantRows <- function(vec, df = intersects_df){
  sapply(1:dim(df)[1],function(x){
    all(df[x,1:2] %in% vec)
  })
}

intersects_df$focal <- relevantRows(SceloporusRanges$species[SceloporusRanges$focal])
intersects_df$Riv_focal <- relevantRows(vec = SceloporusRanges$species[SceloporusRanges$Riv_focal])
intersects_df$definite_overlaps <- intersects_df$GARD & intersects_df$Rivera

#mark the rows in which size classes overlap
intersects_df$samepairs <- NA
intersects_df$samesize <- 

intersects_df$samepairs[relevantRows(
  SceloporusRanges$species[SceloporusRanges$Size_guild == "small"])] <- "small"

intersects_df$samepairs[relevantRows(
  SceloporusRanges$species[SceloporusRanges$Size_guild == "medium"])] <- "medium"

intersects_df$samepairs[relevantRows(
  SceloporusRanges$species[SceloporusRanges$Size_guild == "large"])] <- "large"

intersects_df$samesize <- !is.na(intersects_df$samepairs)

#for any pair of species names, find the corresponding row in intersects_df
pairIndex <- function(vec , df = intersects_df){
  vv <- sapply(1:dim(df)[1], function(x){
    all(vec %in% df[x,1:2])
  })
  if(any(vv)){return(which(vv))}
  else{return(NA)}
  
  
}

saveRDS(intersects_df, file = "Outputs/intersects_df.RDS")

#### Compare our data & Rivera et al's data ####
#which species are in the chromosomally diverse clade and in both datasets?
#trim down the datasets to just the focal species in the Rivera dataset

comparisonsDF <- intersects_df[which(intersects_df$Riv_focal),] 

#where do they disagree?
comparisonsDF$dispute <- comparisonsDF$GARD != comparisonsDF$Rivera
comparisonsDF$dispute <- comparisonsDF$GARD != comparisonsDF$Rivera

comparisonsDF$strictOverlaps <- comparisonsDF$GARD & comparisonsDF$Rivera

comparisonsDF$maxOverlaps <- comparisonsDF$GARD | comparisonsDF$Rivera

saveRDS(comparisonsDF, "Outputs/GARD_Riv_comparisons.RDS")

#Among the focal taxa found in both datasets, where are there disagreements in geographic overlap?
length(which(!comparisonsDF$dispute)) #in agreement on 993 pairs
length(which(comparisonsDF$dispute)) #88 disputed overlaps
length(which(comparisonsDF$GARD[which(comparisonsDF$dispute)])) # GARD dataset estimaes 56 overlaps
length(which(comparisonsDF$Rivera[which(comparisonsDF$dispute)])) #Rivera estimates 32 overlaps
#which taxa are involved in the most disagreements?

problem_species <- sort(table(unlist(comparisonsDF[comparisonsDF$dispute,1:2])), decreasing = TRUE)
# looks like S. undulatus has the most disagreement

#From which dataset do those disagreements spring (which dataset overestimates overlap?)

GARD_problem_species <- sort(table(unlist(
  comparisonsDF[comparisonsDF$dispute & comparisonsDF$GARD,1:2])), decreasing = TRUE)

Rivera_problem_species <- sort(table(unlist(
  comparisonsDF[comparisonsDF$dispute & comparisonsDF$Rivera,1:2])), decreasing = TRUE)



#####Sister species range relationships####
library("phangorn")
library("phytools")

phy <- readRDS(file = "Outputs/SceloporusTree.rds")
intersects_df <- readRDS("Outputs/intersects_df.RDS")
##Find sister species

dd<-lapply(1:phy$Nnode+Ntip(phy),function(n,t)
  Descendants(t,n)[[1]],t=phy)
nodes<-c(1:phy$Nnode+Ntip(phy))[which(sapply(dd,
                                             length)==2)]
sisters<-t(sapply(nodes,function(n,t) 
  t$tip.label[Descendants(t,n)[[1]]],t=phy))

sisters <- gsub("Sceloporus_","",sisters)
sisters<- sisters[-c(8,9,11,12),] 

#find these pairs in the intersects dataframe
sister_int_ind <- sapply(1:dim(sisters)[1], function(x){
  pairIndex(sisters[x,1:2])
})

#make data frame of sister pairs and their intersections
sisters <- as.data.frame(sisters[which(!is.na(sister_int_ind)),])

sisters <- cbind(sisters, intersects_df[sister_int_ind[which(!is.na(sister_int_ind))],
                                        3:5])

# get node depths of sister species pairs

node_depth<-function(phy,node=NULL){ # from Liam Revell (phytools)
  if(is.null(node)) node<-1:phy$Nnode+Ntip(phy)
  h<-max(nodeHeights(phy))
  setNames(h-sapply(node,nodeheight,tree=phy),node)
}


sisters$depths <- node_depth(phy,  nodes)[which(!is.na(sister_int_ind))]

sisters$Karyotype <- as.factor(sisters$Karyotype)
sisters$GARD <- as.factor(sisters$GARD)
sisters$Rivera <- as.factor(sisters$Rivera)

# chi-squared test - can't do it, not enough examples in each category?

library('lsr')
associationTest( ~Karyotype + GARD, sisters)

# Bayesian association test

library('BayesFactor')
crosstab <- xtabs( ~kt_dif + geo_overlap, sisters_testing)
contingencyTableBF( crosstab , sampleType = "jointMulti" ) # BF 1.425053
#A test of association produced a Bayes factor of 1.4:1 in favor of a relationship between karyotype difference and geographic range overlap.



# difference in node depths between overlapping & non-overlapping species
library('ggplot2')
ggplot(sisters_testing)+
  geom_histogram(aes(x= depths, fill = geo_overlap), position="dodge", binwidth=2)+
  geom_vline(xintercept = sisters_testing$depths[which(sisters_testing$kt_dif == "TRUE")])+
  ggtitle("Distribtion of node depths between sister species")


#####analyze species/karyotype overlap####
intersects_df <- readRDS("Outputs/intersects_df.RDS")

intersects_df$Karyotype <- as.factor(intersects_df$Karyotype)
intersects_df$GARD <- as.factor(intersects_df$GARD)
intersects_df$Rivera <- as.factor(intersects_df$Rivera)
intersects_df$definite_overlaps <- as.factor(intersects_df$definite_overlaps)
intersects_df$samesize <- as.factor(intersects_df$samesize)
#Chi sq correlation between Karyotype and Geographic overlap

library('lsr')
chi.GARD.focal.size <- associationTest( ~samesize + GARD, intersects_df[which(intersects_df$focal),])

chi.GARD.full <- associationTest( ~Karyotype + GARD, intersects_df)
chi.GARD.focal <- associationTest( ~Karyotype + GARD, intersects_df[which(intersects_df$focal),])
chi.GARD.riv <- associationTest( ~Karyotype + GARD, intersects_df[which(intersects_df$Riv_focal),])
chi.GARD.focal.med <- associationTest( ~Karyotype + GARD, 
                                       intersects_df[which(intersects_df$focal & intersects_df$samepairs == "medium"),])
chi.GARD.focal.large <- associationTest( ~Karyotype + GARD, 
                                       intersects_df[which(intersects_df$focal & intersects_df$samepairs == "large"),])

chi.Riv.full <-  associationTest( ~Karyotype + Rivera, intersects_df[which(!is.na(intersects_df$Rivera)),])
chi.Riv.focal <- associationTest( ~Karyotype + Rivera, intersects_df[which(intersects_df$Riv_focal),])

def.focal <- associationTest( ~Karyotype + definite_overlaps, intersects_df[which(intersects_df$Riv_focal),])


