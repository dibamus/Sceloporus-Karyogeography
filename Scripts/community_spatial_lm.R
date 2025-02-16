#### COMMUNITY-LEVEL SPATIAL LMS WITH SPATIALREG ####

#Code borrowed from the community_analysis code, updated to run errorSAR models 
community_spatial_lm <- function(rangedf, sizeclass = "kt_div", nresamples = 500){
  library(spatialreg)
  library(spdep)
  
  CRS <- st_crs(rangedf)
  extent <- st_bbox(rangedf)
  
  #Set up grid 
  #r is a raster object used for tallying species and groups
  #grid is an sf dataframe which corresponds to the raster grid r
  #0.5 degree resolution
  r <- raster(extent(matrix( c(extent), nrow=2)), resolution = c(0.5,0.5), 
              crs = CRS$wkt) #"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  grid <- st_as_sf(as(r,'SpatialPolygonsDataFrame'))
  
  spatial_info <- rasterToPoints(r, spatial = TRUE) %>% data.frame()
  
  colnames(spatial_info) <- c("lon","lat")
  
  #Check for range intersections between Sceloporus species and grid squares
  #spit out sparse matrix
  speciesgrid <- st_intersects(grid,rangedf, sparse = TRUE)
  
  #species count values for pixels in raster
  spatial_info$sp_values <- lengths(speciesgrid) #
  
  #create neighbor list for spatialreg ----
  nblist <- cell2nb(nrow = nrow(r), ncol = ncol(r), type="queen", torus=FALSE, legacy=FALSE, x=NULL)
  
  set.ZeroPolicyOption(TRUE)
  
  nb <- nb2listw(nblist, zero.policy = TRUE)
  
  
  #Translate species counts to karyotype counts ----
  species_indices <- unlist(speciesgrid) #x is a list of integers. each element corresponds to an instance of a sceloporus species in a grid square. each unique integer corresponds to a species
  
  #the next steps look up which species are in which karyotype groups and then substitute karyotype group for species in x
  
  from <- c(1:length(rangedf$chromosome_group))
  
  to <- rangedf$chromosome_group #this is the vector that needs to be resampled for frequentist comparisons
  map <- setNames(to, from) #now we can look up karyotype group based on the species integers in x
  
  species_indices[] <- map[species_indices]
  
  karyotypegrid <- relist(species_indices, skeleton=speciesgrid)
  
  ##=====##
  #Karyotype count values for pixels in raster
  spatial_info$kt_values <- sapply(karyotypegrid, function(x) length(unique(x))) #count the number of unique karyotypes in each grid square
  
  # #Remove grid squares with 0 or 1 species in them; these will necessarily have 0 or 1 karyotypes
  # replace_1_0 <- function(x){
  #   x[which(x$sp_values < 2),3:4]<- NA
  #   return(x)
  # }
  # spatial_info <- replace_1_0(spatial_info)
  
  set.ZeroPolicyOption(TRUE)
  errorSAR <- errorsarlm(kt_values ~ sp_values, data = spatial_info, listw = nb, 
                         na.action = na.exclude, zero.policy = T) 
  
  errorSAR$resampled_karyotypes <- FALSE

  
  ##RESAMPLING
  #quickly compare actual data correlation between species and karyotype to 
  #simulated data where karyotypes are distributed randomly to species
  
  gp <- rangedf$chromosome_group
  
  x <- unlist(speciesgrid)
  
  from <- c(1:length(gp))
  
  # Resample karyotypes ----
  kt.res <-function(gp = rangedf$chromosome_group, 
                    x = unlist(speciesgrid), 
                    from  = c(1:length(gp))){
    to <- sample(gp) #this is the vector that needs to be re-sampled for frequentist comparisons
    map <- setNames(to, from)
    
    x[] <- map[x]
    
    ktres <- relist(x, skeleton=speciesgrid)
    ktval <- sapply(ktres, function(x) length(unique(x)))
    return(ktval)
  }
  #correlate resampled grids ----
  correlate.resamples <- function(df = spatial_info, listw = nb, kt_resampler = kt.res){
    df$kt_values <- kt_resampler()
    model <- spatialreg::errorsarlm(
      kt_values ~ sp_values, data = df, #making 0 and 1 values NA
      listw = listw, na.action = na.exclude, zero.policy = T)
    model$resampled_karyotypes <- TRUE
    return(model)
  }
  library("parallel")
  library("foreach")
  cl <- parallel::makeCluster(3)
  doParallel::registerDoParallel(cl)
  
  sim_list <- foreach(i = 1:nresamples) %dopar% correlate.resamples(df = spatial_info, 
                                                           listw = nb, 
                                                           kt_resampler = kt.res)
  
  parallel::stopCluster(cl)
  
  saveRDS(sim_list, paste0("Results/spatial/",sizeclass,"_spatial_lm_sims.rds"))
  saveRDS(errorSAR, paste0("Results/spatial/",sizeclass,"_spatial_lm_real.rds"))
}