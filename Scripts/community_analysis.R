#Karyo-Community
#Within communities, how are species diversity and karyotypic diversity related?
#Outputs
# 
# a raster stack with 3 layers
# 1 - species diversity
# 2 - karyotype diversity
# 3 - species diversity - karyotype diversity
#
# 1 data frame of species diversity & karyotype diversity in each grid cell
#
# 1 set of 1000 resampled karyotype-species mappings
# 1 plot of the distribution of possible species/karyptype correlations 
# versus the observed species/karyotype diversity

#### COMMUNITY-LEVEL ANALYSIS ####
#Copy geographic data from Sceloporus Ranges

community_analysis <- function(rangedf, sizeclass = "kt_div", nresamples = 10000){

CRS <- st_crs(rangedf)
extent <- st_bbox(rangedf)


#Set up grid 
#r is a raster object used for tallying species and groups
#grid is an sf dataframe which corresponds to the raster grid r
#0.5 degree resolution
r <- raster(extent(matrix( c(extent), nrow=2)), resolution = c(0.5,0.5), 
            crs = CRS$wkt) #"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
grid <- st_as_sf(as(r,'SpatialPolygonsDataFrame'))


#Check for range intersections between Sceloporus species and grid squares
#spit out sparse matrix
speciesgrid <- st_intersects(grid,rangedf, sparse = TRUE)


#Translate species counts to karyotype counts
x <- unlist(speciesgrid) #x is a list of integers. each element corresponds to an instance of a sceloporus species in a grid square. each unique integer corresponds to a species

#the next steps look up which species are in which karyotype groups and them substitute karyotype group for species in x

from <- c(1:length(rangedf$chromosome_group))

to <- rangedf$chromosome_group #this is the vector that needs to be resampled for frequentist comparisons
map <- setNames(to, from) #now we can look up karyotype group based on the species integers in x

x[] <- map[x]

karyotypegrid <- relist(x, skeleton=speciesgrid)

##=====##
#Karyotype count values for pixels in raster
kt_values <- sapply(karyotypegrid, function(x) length(unique(x))) #count the number of unique karyotypes in each grid square

kt_commons <- sapply(karyotypegrid, function(x) sort(table(x),decreasing=TRUE)) # not currently used for anything: makes a list of the most common karyotypes (descending order of commonness) in each grid square


#species count values for pixels in raster
sp_values <- lengths(speciesgrid)

# Generate raster layers for sceloporus species and karyotype diversity
#species diversity raster
r_sp <- r #copy raster grid
r_sp[] <- sp_values #value of each grid square corresponds to the # of species in that square

sp_table <- table(sp_values) #check out how many squares have x number of species

#karyotype diversity raster
r_kt <- r
r_kt[] <- kt_values

kt_table <- table(kt_values)

#difference between r_sp and r_kt
r_diff <- r_sp - r_kt 

excess_kt_sp <- r_diff/r_sp


##VECTOR COMPARISON
# since all pixels with 0 and 1 species necessarily contain 0 and 1 karyotype groups, they cannot be included in an analysis of correlation between species diversity and karyotype diversity
#therefore, correlating species and karyotype vectors needs a custom function

KaryoCorr <- function(sp,kt, range) { #0.7516846 for 1 degree squares, 0.7010471 for 0.5 degree squares
  #range argument should be set to c(2:11) for regular results
  #setting other variables for range argument allows for finer understanding of the species-karyotype relationship at different species diversity levels.
  if (length(sp) != length(kt)){
    print("vectors of unequal lengths")
    return()
  }
  else {
    x <- which(sp %in% range)
    c(cor(sp[x],kt[x]), # Pearson's R
      lm(kt[x] ~ sp[x])$coefficients[2]) # slope
  }
}

range <- c(2:11)# set focal diversity range

kt.cor <- KaryoCorr(sp_values,kt_values,range)

##RESAMPLING
#quickly compare actual data correlation between species and karyotype to 
#simulated data where karyotypes are distributed randomly to species

gp <- rangedf$chromosome_group

x <- unlist(speciesgrid)

from <- c(1:length(gp))

kt.res <-function(){
  to <- sample(gp) #this is the vector that needs to be re-sampled for frequentist comparisons
  map <- setNames(to, from)
  
  x[] <- map[x]
  
  ktres <- relist(x, skeleton=speciesgrid)
  ktval <- sapply(ktres, function(x) length(unique(x)))
  return(ktval)
}

correlate.resamples <- function(){
  rkt <- kt.res()
  return(KaryoCorr(sp_values,rkt,range))
}

resampled_ktvals <- replicate(nresamples, correlate.resamples()) %>% t()
colnames(resampled_ktvals) <- c("R","slope")

saveRDS(resampled_ktvals, paste0("Results/",sizeclass, "/", sizeclass, "resampledktvals.rds"))

Obs.cor <- KaryoCorr(sp_values,kt_values, range)

kt.cor.p <- length(resampled_ktvals[,1][which(resampled_ktvals[,1] < Obs.cor[1])])/nresamples #(this is a p value) (p= 0.6494 with "resampledktvals-05x05.R")
kt.slope.p <- length(resampled_ktvals[,2][which(resampled_ktvals[,2] < Obs.cor[2])])/nresamples #(this is a p value) (p= 0.6494 with "resampledktvals-05x05.R")

#save this correlation result
saveRDS(list(cor = Obs.cor[1], 
             slope = Obs.cor[2], 
             cor.p = kt.cor.p,
             slope.p = kt.slope.p),paste0("Results/",sizeclass, "/", sizeclass, "_observed_correlation.rds"))

# save the rasters
rStack <- stack(r_sp,r_kt,r_diff)
names(rStack) <- c('species','karyotype','difference')

saveRDS(rStack, paste0("Results/",sizeclass, "/", sizeclass, "_rasters.rds"))

#save the observed sp & kt values
gridDF <- data.frame(kt_values,sp_values)
saveRDS(gridDF,paste0("Results/",sizeclass, "/", sizeclass, "_observedDF.rds"))

}