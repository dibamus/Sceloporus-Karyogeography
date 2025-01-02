#### PLOTS ####

library("ggplot2")
library("sf")
library("raster")
library("terra")
library('tidyterra')
library("rnaturalearth")
fullset <- readRDS("Results/kt_div/kt_div_rasters.rds") %>% terra::rast()
small <- readRDS("Results/small/small_rasters.rds") %>% terra::rast()
medium <- readRDS("Results/medium/medium_rasters.rds") %>% terra::rast()
large <- readRDS("Results/large/large_rasters.rds") %>% terra::rast()

world <- st_as_sf(ne_countries(scale = "medium", returnclass = "sf"))

scel_aes <- c(scale_x_continuous(limits = ext(fullset)[c(1,2)], expand = c(0, 0)),
              scale_y_continuous(limits = ext(fullset)[c(3,4)], expand = c(0, 0)))

scel_theme <- theme(axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    axis.ticks=element_blank(),
                    strip.text.x = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    legend.position = "none")

left_theme <- theme(axis.title.y= element_text(size = 12,
                                               angle = 90))
right_theme <- theme(strip.text.x = element_text(),
                    strip.text.y = element_blank(),
                    axis.title.y= element_blank(),
                    legend.position = "right",
                    legend.position.inside = c(1,1),
                    legend.justification = c(1,1),
                    legend.title = element_blank(),
                    legend.direction="vertical",
                    legend.text=element_text(size=10),
                    legend.text.position = "left")
                    #legend.background = element_rect(fill = "#ffffff"))

fsmax <- as.character(minmax(fullset$species)[2])
fsmax_kt <- as.character(minmax(fullset$karyotype)[2])
fsdiff <- as.character(minmax(fullset$species - fullset$karyotype)[2])

spgrad <- scale_fill_gradient(low ="#7a7a7a00",high ="#10ccdaff", na.value="#00000000", 
                              breaks = c(1, minmax(fullset$species)[2]), 
                              limits = c(1, minmax(fullset$species)[2]),
                              labels = c("0",fsmax),
                              name = "")

#Sceloporus community diversity maps
#full set
fullset_plot_sp <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = fullset$species) +
  spgrad +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", fsmax)),
            aes(x = x, y = y, label = text), size = 3) +
  
  scel_aes +
  scel_theme +
  left_theme +
  ylab("species diversity") +
  
  ggtitle('all species')
  

fullset_plot_kt <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = fullset$karyotype) +
  ylab("karyotype diversity") +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ffdd66ff", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$karyotype)[2]), 
                      limits = c(1, minmax(fullset$karyotype)[2]), 
                      labels = c("0",fsmax_kt),
                      name = "") +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", fsmax_kt)),
            aes(x = x, y = y, label = text), size = 3) +
  
  scel_aes +
  scel_theme +
  left_theme

fullset_plot_diff <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = fullset$difference) +
  ylab("difference") +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ff7777ff", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$difference)[2]), 
                      limits = c(1, minmax(fullset$difference)[2]), 
                      labels = c("0",fsdiff),
                      name = "") +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           fsdiff)),
            aes(x = x, y = y, label = text), size = 3) +
  
  scel_aes +
  scel_theme +
  left_theme


#plot smalls

small_plot_sp <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = small$species) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#10ccdaff", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$species)[2]), 
                      limits = c(1, minmax(fullset$species)[2]),
                      labels = c("0",fsmax)) +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(small$species)[2])),
            aes(x = x, y = y, label = text), size = 3) +
  scel_aes +
  scel_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())+
  ggtitle('small species')

small_plot_kt <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = small$karyotype) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ffdd6688", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$karyotype)[2]), 
                      limits = c(1, minmax(fullset$karyotype)[2]),
                      labels = c("0",fsmax)) +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(small$karyotype)[2])),
            aes(x = x, y = y, label = text), size = 3) +
  scel_aes +
  scel_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())

small_plot_diff <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = small$difference) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ff777725", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$difference)[2]), 
                      limits = c(1, minmax(fullset$difference)[2]),
                      labels = c("0",fsmax)) +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(small$difference)[2])),
            aes(x = x, y = y, label = text), size = 3) +

  scel_aes +
  scel_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())

#plot mediums

medium_plot_sp <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = medium$species) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#10ccda8c", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$species)[2]),  
                      limits = c(1, minmax(fullset$species)[2]),
                      labels = c("0",fsmax)) +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(medium$species)[2])),
            aes(x = x, y = y, label = text), size = 3) +

  scel_aes +
  scel_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())+
  ggtitle('mid-sized species')

medium_plot_kt <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = medium$karyotype) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ffdd66aa", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$karyotype)[2]),  
                      limits = c(1, minmax(fullset$karyotype)[2]),
                      labels = c("0",fsmax)) +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(medium$karyotype)[2])),
            aes(x = x, y = y, label = text), size = 3) +

  scel_aes +
  scel_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())

medium_plot_diff <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = medium$difference) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ff77776e", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$difference)[2]),  
                      limits = c(1, minmax(fullset$difference)[2]),
                      labels = c("0",fsmax)) +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(medium$difference)[2])),
            aes(x = x, y = y, label = text), size = 3) +

  scel_aes +
  scel_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())

#larges

large_plot_sp <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = large$species) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#10ccdaaa", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$species)[2]), 
                      limits = c(1, minmax(fullset$species)[2]), 
                      labels = c("0",fsmax)) +
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(large$species)[2])),
            aes(x = x, y = y, label = text), size = 3) +

  scel_aes +
  scel_theme +
  right_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())+
  ggtitle('large species')

large_plot_kt <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = large$karyotype) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ffdd66a2", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$karyotype)[2]), 
                      limits = c(1, minmax(fullset$karyotype)[2]), 
                      labels = c("0",fsmax_kt)) +
  
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(large$karyotype)[2])),
            aes(x = x, y = y, label = text), size = 3) +

  scel_aes +
  scel_theme +
  right_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())

large_plot_diff <- ggplot() +
  geom_sf(data = world, fill = "white")+ theme_bw()+
  geom_spatraster(data = large$difference) +
  scale_fill_gradient(low ="#7a7a7a00",high ="#ff777792", na.value="#00000000", 
                      breaks = c(1, minmax(fullset$species-fullset$karyotype)[2]), 
                      limits = c(1, minmax(fullset$species-fullset$karyotype)[2]), 
                      labels = c("0",as.character(minmax(fullset$species-fullset$karyotype)[2]))) +
  
  geom_text(data = data.frame(x = -113, y = 12, 
                              text = paste("max = ", 
                                           minmax(large$difference)[2])),
            aes(x = x, y = y, label = text), size = 3) +
  scel_aes +
  scel_theme +
  right_theme +
  theme(strip.text.x = element_text(),
        strip.text.y = element_text())

#correlation data

cor_aes <- c(scale_x_continuous(limits = c(0,1), expand = c(0, 0)),
             scale_y_continuous(limits = c(0,6), expand = c(0,0),)
             )
cor_theme <- theme_bw() + theme(panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.title.x = element_text("slope"),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   aspect.ratio= 0.5)

#follow this example
kt_df <- readRDS("Results/kt_div/kt_divresampledktvals.rds") %>% as.data.frame()
kt.cor <- readRDS("Results/kt_div/kt_div_observed_correlation.rds")

kt_dens <- ggplot(kt_df, aes(x=slope)) + 
  geom_density(fill = "#ccc") +
  geom_vline(xintercept = kt.cor$slope, color = "#10ccda", 
             linetype = "dotted",
             linewidth = 1) + cor_aes + cor_theme + 
  theme(axis.title.y = element_text("density", angle = 90))


small_df <- readRDS("Results/small/smallresampledktvals.rds") %>% as.data.frame()
small.cor <- readRDS("Results/small/small_observed_correlation.rds")
small_dens <- ggplot(small_df, aes(x=slope)) + 
  geom_density(fill = "#ccc") +
  geom_vline(xintercept = small.cor$slope, color = "#10ccda", 
linetype = "dotted",
linewidth = 1) + cor_aes + cor_theme

medium_df <- readRDS("Results/medium/mediumresampledktvals.rds") %>% as.data.frame()
medium.cor <- readRDS("Results/medium/medium_observed_correlation.rds")
medium_dens <- ggplot(medium_df, aes(x=slope)) + 
  geom_density(fill = "#ccc") +
  geom_vline(xintercept = medium.cor$slope, color = "#10ccda", 
             linetype = "dotted",
             linewidth = 1) + cor_aes + cor_theme

large_df <- readRDS("Results/large/largeresampledktvals.rds") %>% as.data.frame()
large.cor <- readRDS("Results/large/large_observed_correlation.rds")
large_dens <- ggplot(large_df, aes(x=slope)) + 
  geom_density(fill = "#ccc") +
  geom_vline(xintercept = large.cor$slope, color = "#10ccda", 
             linetype = "dotted",
             linewidth = 1) + cor_aes + cor_theme


library('patchwork')

big_grid <- (fullset_plot_sp / fullset_plot_kt / fullset_plot_diff /kt_dens) |
  (small_plot_sp / small_plot_kt / small_plot_diff /small_dens) |
  (medium_plot_sp / medium_plot_kt / medium_plot_diff /medium_dens) |
  (large_plot_sp / large_plot_kt / large_plot_diff /large_dens)
dev.off()
gc()
big_grid

ggsave("Plots/Maps&Correlations.pdf",big_grid, width = 8, height = 8, units = "in")
ggsave("Plots/Maps&Correlations.png",big_grid, width = 8, height = 8, units = "in")
