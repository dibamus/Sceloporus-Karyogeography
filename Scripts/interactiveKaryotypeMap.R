#Sceloporus KaryoPolygons for Figure 1
library('sf')
library('sp')
library('dplyr')
library('ggplot2')
library('leaflet')
library('colorspace')
library('htmlwidgets')
ranges <- readRDS("Data/Sceloporus_data_frame.RDS")

kranges <- ranges %>%
  group_by(chromosome_group) %>% 
  reframe(geometry = st_union(geometry))
area <- st_area(st_make_valid(kranges$geometry))

kranges <- kranges[order(area),]

factpal <- colorFactor(qualitative_hcl(17), kranges$chromosome_group)

map <- leaflet(kranges) %>%
  addTiles() %>%
  addPolygons(data = kranges$geometry, 
              color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.2,
              fillColor = factpal(kranges$chromosome_group),
              label = kranges$chromosome_group,
              highlightOptions = highlightOptions(color = "white", weight = 2,
                                                  bringToFront = TRUE))

saveWidget(map, file = 'Plots/interactivemap.html')
