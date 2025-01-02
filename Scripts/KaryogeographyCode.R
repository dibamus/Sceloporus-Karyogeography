#### Sceloporus Karyogeography
####create sceloporus geographic dataset####
setwd("C:/Users/Isaac/Syncthing-Docs/Sceloporus Karyogeography")

library("sf")
library("dplyr")
#Load in the Sceloporus from the GARD dataset (original dataset not included)
load("Data/SceloporusRanges.sf")
#Load in a CSV of species-specific data
species.data <- read.csv("Data/Chromosomes.csv", fileEncoding = 'UTF-8-BOM')

row.names(species.data) <- species.data$Binomial

species.data$tree_name <- gsub(" ","_",species.data$Binomial)

species.data$species <- gsub("Sceloporus_","",species.data$tree_name)

species.data$focal <- species.data$chromosomally_diverse_clade &
  species.data$chromosome_group != "" &
  species.data$insular_endemic == FALSE
species.data$Riv_focal <- species.data$focal & species.data$in_Rivera_et_al

SceloporusRanges <- left_join(SceloporusRanges, species.data, by = "Binomial")



saveRDS(SceloporusRanges, "Outputs/Sceloporus_data_frame.RDS")

#### Community Analysis####
library("sf")
library("dplyr")
library("raster")
source("community_analysis.R")
SceloporusRanges <- readRDS("Outputs/Sceloporus_data_frame.RDS")
#remove the karyotypically non-diverse clade & the island endemics

SceloporusRanges_focal <- subset(SceloporusRanges, focal)

community_analysis(SceloporusRanges,"kt_div")

#For SMALL
small.data <- subset(SceloporusRanges_focal, Size_guild == "small")
community_analysis(small.data,"small")

#for MEDIUM
medium.data <- subset(SceloporusRanges_focal, Size_guild == "medium")
community_analysis(medium.data,"medium")

#for LARGE
large.data <- subset(SceloporusRanges_focal, Size_guild == "large")
community_analysis(large.data,"large")

# #### Spatial Correlation ####
# library("SpatialPack")
# spts <- rasterToPoints(r_sp, spatial = TRUE)
# SPdf <- data.frame(spts@coords[,1],
#                    spts@coords[,2],
#                    sp_values,
#                    kt_values)
# 
# SPdf<- SPdf[which(SPdf$sp_values >1),]
# 
# SpatialCorr <- cor.spatial(x = SPdf$sp_values,
#                            y = SPdf$kt_values, 
#                            coords = SPdf[,c(3,4)])
# 
# SpatialtTest <- modified.ttest(x = SPdf$sp_values,
#                            y = SPdf$kt_values, 
#                            coords = SPdf[,c(3,4)])

#### Size Guilds VS Karyotypes ####
#Are Sceloporus size guilds more karyotypically diverse than expected by chance?
species.data <- readRDS("Data/Sceloporus_data_frame.rds")
species.data <- filter(species.data, focal == TRUE)

guild_n <- table(species.data$Size_guild)

resample_sizes <- function(df, guild, n = 10000){
  table <- filter(df, Size_guild == guild)$chromosome_group %>% table
  total_sp <- sum(table)
  total_chr <- length(table)
  sims <- rep(NA, times = n)
  for (i in 1:n){
    sims[i] <-length(unique(df[sample(dim(df)[1],total_sp),]$chromosome_group))
  }
  
  return(list(sims = sims, observed = total_chr))
}

sm <- resample_sizes(species.data, guild = "small")
sm.p <- length(which(sm$sims <= sm$observed))/10000 #0.71
med <- resample_sizes(species.data, guild = "medium")
med.p <- length(which(med$sims <= med$observed))/10000 #0.17
lg <- resample_sizes(species.data, guild = "large")
lg.p <- length(which(lg$sims <= lg$observed))/10000 #0.18

#### PHYLOGENY ####
library("picante")
#### Adjust phylogeny ####
phy2 <- read.tree(file = "Data/Scel_tree.txt")

#Sceloporus tips to remove:
rm <- c("Sceloporus_occidentalis_uwbm6281",
        "Sceloporus_formosus_rvt76",
        "Sceloporus_grammicus_microlepidotus",
        "Sceloporus_formosus_scitulus",
        "Sceloporus_megalepidurus_pictus",
        "Sceloporus_torquatus_uogv2526",
        "Sceloporus_insignis_anmo1130",
        "Sceloporus_aureolus_rvt54",
        "Sceloporus_orcutti_uwbm7654",
        "Sceloporus_magister_uniformis_mvz162077",
        "Sceloporus_zosteromus_ADG74",
        "Sceloporus_samcolemani_rwb06263",
        "Sceloporus_scalaris_rwb06247",
        "Sceloporus_subniger_rwb0645",
        "Sceloporus_pyrocephalus_utar53473",
        "Sceloporus_siniferus_uwbm6653"
)
#Sceloporus tips to rename:
#current names:
oldnames <- c("Sceloporus_occidentalis_mvz245697",
"Sceloporus_acanthinus_anmo1932",
"Sceloporus_taeniocnemis_mvz4213",
"Sceloporus_formosus_anmo1248",
"Sceloporus_torquatus_uwbm6600",
"Sceloporus_serrifer_utar39870",
"Sceloporus_aureolus_jac22409",
"Sceloporus_orcutti_rwm798",
"Sceloporus_magister_bimaculosus_dgm924",
"Sceloporus_magister_uniformis_dgm474",
"Sceloporus_magister_cephaloflavus_uwbm7395",
"Sceloporus_magister_magister_rom14488",
"Sceloporus_zosteromus_ADG49",
"Sceloporus_samcolemani_jjw698",
"Sceloporus_scalaris_uwbm6589",
"Sceloporus_subniger_rwb0686",
"Sceloporus_pyrocephalus_unknown",
"Sceloporus_siniferus_mvz236299",
"Sceloporus_tristicus",
"Sceloporus_cowelsi"
)
#new names:
newnames<-c("Sceloporus_occidentalis",
"Sceloporus_acanthinus",
"Sceloporus_taeniocnemis",
"Sceloporus_formosus",
"Sceloporus_torquatus",
"Sceloporus_serrifer",
"Sceloporus_aureolus",
"Sceloporus_orcutti",
"Sceloporus_bimaculosus",
"Sceloporus_uniformis",
"Sceloporus_cephaloflavus",
"Sceloporus_magister",
"Sceloporus_zosteromus",
"Sceloporus_samcolemani",
"Sceloporus_scalaris",
"Sceloporus_subniger",
"Sceloporus_pyrocephalus",
"Sceloporus_siniferus",
"Sceloporus_tristichus",
"Sceloporus_cowlesi"
)


phy2 <- drop.tip(phy2,rm)

matches<-match(oldnames,phy2$tip.label)

phy2$tip.label[matches] <- newnames

nonSceltips <- grep("Sceloporus",phy2$tip.label,invert = TRUE)

phy2 <- drop.tip(phy2,nonSceltips)

saveRDS(phy2,file = "Outputs/SceloporusTree.rds")

#### Plot Simmap ####
library("phytools")
library("tidyverse")
species.data <- readRDS("Outputs/Sceloporus_data_frame.rds")
species.data <- species.data[which(species.data$chromosome_group !=""),]
phy <- readRDS(file = "Outputs/SceloporusTree.rds")

species.data$treeName <- gsub(" ","_",species.data$Binomial)

phy <- drop.tip(phy, which(!(phy$tip.label %in% species.data$treeName)))

species.data$treeName[which(!(species.data$treeName %in% phy$tip.label))]

chromosomegroup <- species.data$chromosome_group[species.data$treeName %in% phy$tip.label] %>% as.factor()
names(chromosomegroup)<- species.data$treeName[species.data$treeName %in% phy$tip.label]


ScSimMap <- make.simmap(phy,chromosomegroup,model = "ER", nsim = 100)#

cols <- setNames(c("darkslategray3","gray60","gold","lightskyblue","darkolivegreen3","tomato","orange1","orangered3","orchid","mediumpurple","darkorchid","black","royalblue4","rosybrown"),
                 levels(chromosomegroup))
plot(ScSimMap[[1]],cols,fsize=0.8,lwd=4,ftype="i")
add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE,fsize=0.9)

plot(summary(ScSimMap),colors=cols,fsize=0.9,cex=c(0.5,0.2))
add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE,fsize=0.9)



#3 pairs of sister species have different chromosome arrangements
diffSis <- c("Sceloporus_zosteromus","Sceloporus_lineatulus",
             "Sceloporus_palaciosi","Sceloporus_anahuacus",
             "Sceloporus_maculosus","Sceloporus_gadoviae") %>%
  matrix(nrow =3,ncol =2,byrow = T)

which(sisters %in% diffSis)

