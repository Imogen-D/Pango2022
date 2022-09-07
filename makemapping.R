#!/usr/bin/env Rscript

#module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R



# call libraries 
library(raster)
library(rgdal)
library(ggplot2)
library(reshape2)
library(plyr)
library(scatterpie)
library(dplyr)
library(rgdal)
library(raster) 
library(rworldmap) 
library(ggrepel)
library(RColorBrewer)
library(raster)
library(magick)
library(cowplot)

setwd('/work/idumville/pangolins/RAD/data/geodata')

##########################################
####### Plotting with shape files  ########## 
##########################################

                                      
ne_rivers <- readOGR('../../data/geodata/ne_10m_rivers_lake_centerlines.shp')
ne_count <- readOGR('../../data/geodata/ne_10m_admin_0_countries.shp')
ne_roads <- readOGR('../../data/geodata/ne_10m_roads.shp')
ne_rail <- readOGR('../../data/geodata/ne_10m_railroads.shp')
ne_urban <- readOGR('../../data/geodata/ne_10m_urban_areas.shp')
ne_lakes <- readOGR('../../data/geodata/ne_10m_lakes.shp')
#ne_ocean <- readOGR('../../data/geodata/ne_10m_coastline.shp') #ocean.shp') #bathymetry_L_0.shp')
#forest <- readOGR('../../data/geodata/ifl_intact_forest_landscapes_v2018.shp')
ne_labels <- readOGR('../../data/geodata/wld_poi_admin0_unmap_2019.shp')#ne_10m_geography_regions_polys.shp')
pango<-readOGR('../../data/geodata/data_0.shp')
 
 
lon <- c(8.000, 14.2194) #-10.5318, 29.53505) #
lat <- c(-0.271891, 7.25) #-3.7346, 10.7618) #
domain <- c(8.000, 14.2194, -0.271891, 7.25) #-10.5318, 29.53505, -3.7346, 10.7618) #

smallpoly <- data.frame(lat =c(8.000, 14.2194, 14.2194, 8.000), long= c(-0.271891, -0.271891, 7.25, 7.25))

bigpoly <- data.frame(lat =c(-10.5318, 12, 12, -10.5318), long= c(-3.7346, -3.7346, 7,7))

                                      
for.earth <- raster('../../data/geodata/Crowther_Nature_Files_Revision_01_WGS84_GeoTiff/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif')

for.crop <- crop(for.earth, extent(domain)) 
for.crop <- aggregate(for.crop, fact=2) #10 for big map, 2 for small

for.table <- as.data.frame(for.crop, xy=T)
  for.table <- na.omit(for.table)
  

river.subset <- crop(ne_rivers, extent(domain))
count.subset <- crop(ne_count, extent(domain))
road.subset <- crop(ne_roads, extent(domain))
rail.subset <- crop(ne_rail, extent(domain))
urban.subset <- crop(ne_urban, extent(domain))
lakes.subset <- crop(ne_lakes, extent(domain))
pango.subset <- crop(pango, extent(domain))
#ocean.subset <- crop(ne_ocean, extent(domain)) 
label.subset <- as.data.frame(crop(ne_labels, extent(domain)))


label.subset <- label.subset[-c(which(label.subset$Label %in% c("Sao Tome and Principe", "Abyei"))),] #sao tome and principe


                                      
p0 <- ggplot() + 
  geom_polygon(data=smallpoly, aes(x=lat, y=long), fill = "cadetblue3") + 
  geom_polygon(data=pango.subset, aes(x=long, y=lat, group=group), colour = 'NA', fill = "antiquewhite", linetype = 3, size = 0.3) +
  geom_tile(data=for.table, aes(x=x, y=y, fill=Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff), show.legend = FALSE) +
  #adding boyndary pango line
  geom_polygon(data=pango.subset, aes(x=long, y=lat, group=group), colour = 'red', fill = NA, linetype = 3, size = 0.3) +
  geom_polygon(data=urban.subset, aes(x = long, y = lat, group = group), fill = "bisque3", colour = "bisque3") +
  geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'deepskyblue2', size=0.3) +
  geom_polygon(data=lakes.subset, aes(x = long, y = lat, group = group), fill = "deepskyblue2", colour = "deepskyblue2") +
  geom_path(data=count.subset, aes(x = long, y = lat, group = group), color = 'black', size=0.3) +
  
 #draw_image("/work/idumville/pangolins/RAD/data/geodata/Manistricuspis-removebg-preview.png", x = 11, y = 5, width = 3, height = 3) + #x = -10, y = -7, width =10, height =10) + #
  draw_image("/work/idumville/pangolins/RAD/data/geodata/plane-removebg-preview.png", x = 8, y = 1, width = 1, height = 1) +  #x = 3, y = -2, width = 2, height =2) + #
  geom_path(data=road.subset, aes(x = long, y = lat, group = group), color = 'azure4',  size=0.3) +
  geom_path(data=rail.subset, aes(x = long, y = lat, group = group), color = 'azure3', size=0.3) +
  scale_fill_gradientn(colors=c(NA,  "#A0D600FF", "#63C600FF","#63C600FF", "#63C600FF")) + #"#00A600FF", "#2DB600FF", only 3 NAs for small, 6 for big + #old = terrain.colors(10) for big
  geom_text(data=label.subset, aes(x = coords.x1, y = coords.x2, label = Label),colour = "black", size=2) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0))  +
  xlab('') + ylab('') + coord_fixed() + theme_bw() #+ xlim(lon) + ylim(lat)

pdf('../../data/geodata/testsmallmapNOPANGONOPLANE.pdf')
plot(p0)
dev.off()

saveRDS(p0, file="../../data/geodata/smallmap.rdata")


