#!/usr/bin/env Rscript

#module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R

library(raster) 
library(sf)
library(rgdal)
library(gdalUtils)

setwd('/work/idumville/pangolins/RAD/data/geodata')

##Cropping
Cam <- extent(c(8.000, 14.2194), c(-0.271891, 7.25))

#for1 <- crop(raster("/work/idumville/pangolins/RAD/data/geodata/Hansen/Hansen_GFC-2021-v1.9_treecover2000_00N_000E.tif"), Cam) #aggregate(,2)
#for2 <- crop(raster("/work/idumville/pangolins/RAD/data/geodata/Hansen/Hansen_GFC-2021-v1.9_treecover2000_00N_010E.tif"), Cam)
#for4 <- crop(raster("/work/idumville/pangolins/RAD/data/geodata/Hansen/Hansen_GFC-2021-v1.9_treecover2000_10N_000E.tif"), Cam)
#for5 <- crop(raster("/work/idumville/pangolins/RAD/data/geodata/Hansen/Hansen_GFC-2021-v1.9_treecover2000_10N_010E.tif"), Cam)

#merge(for1, for2, for4, for5, filename="/work/idumville/pangolins/RAD/data/geodata/Hansen/allforest.tif")



hfp <- raster("./hfp-africa-geo-grid/hfp_Africa_grid/hfp_africa/w001001.adf")
hfp <- crop(hfp, Cam)
lulc <- raster("./lulc-human-modification-terrestrial-systems_geographic.tif")
lulc <- crop(lulc, Cam)
forest <- raster('./Crowther_Nature_Files_Revision_01_WGS84_GeoTiff/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif')
forest <- crop(forest, Cam)
forest2 <- raster("/work/idumville/pangolins/RAD/data/geodata/Hansen/allforest.tif")
forest2 <- crop(forest2, Cam)


#convert points to spatial
#pred <-  read.delim("../../results/locator/indvidualcentroids.txt")
#pred <- read.delim("/work/idumville/pangolins/RAD/results/locator/seeding.concat.allWCA/seedinglocator00/seeding00error_centroids.txt")
pred <- read.delim("/work/idumville/pangolins/RAD/results/locator/seeding.concat.allWCA/seedinglocator90/seeding90error_centroids.txt")

pts <- SpatialPoints(pred[,2:3],  proj4string = CRS(proj4string(forest)))

## compare values

pred$forest <- extract(forest, pred[,2:3]) #dataframe of population and information
pred$lulc <- extract(lulc, pred[,2:3])
pred$hfp <- extract(hfp, pred[,2:3])

pred$forest2 <- extract(forest2, pred[,2:3]) #dataframe of population and information


#intersection of polygons and coordinates
#input shape file for BIOMEs
tt <- read_sf('ltw-africa-geo-grid/ltw_africa.shp')

pnts_sf <- sf::st_as_sf(pred, coords = c("kd_x", "kd_y"), crs = 4326)
pnts_trans <- st_transform(pnts_sf, "EPSG:32632") 
tt_trans <- st_transform(tt, "EPSG:32632")  
res <- sf::st_join(pnts_trans, tt_trans)
print(res)
resdf <- as.data.frame(res)
biomes <- resdf[c(which(resdf$BIOME != "NA")),] #whichparks are they in

library(rgeos)
spts = SpatialPoints(SpatialPoints(pred[,2:3]))
pnts_trans <- as(pnts_trans, "Spatial")
tt_trans <- as(tt_trans, "Spatial")

BIOMEdist <- apply(gDistance(pnts_trans, tt_trans,byid=TRUE),2,min)

pred$biomedist <- BIOMEdist

#input shape file for national parks
unzip("/work/idumville/pangolins/RAD/data/geodata/WDPA_WDOECM_Jun2022_Public_AF_shp_0.zip")
prot1 <- crop(readOGR("/work/idumville/pangolins/RAD/data/geodata/WDPA_WDOECM_Jun2022_Public_AF_shp-polygons.shp"), Cam)

unzip("/work/idumville/pangolins/RAD/data/geodata/WDPA_WDOECM_Jun2022_Public_AF_shp_1.zip")
prot2 <- crop(readOGR("/work/idumville/pangolins/RAD/data/geodata/WDPA_WDOECM_Jun2022_Public_AF_shp-polygons.shp"), Cam)


unzip("/work/idumville/pangolins/RAD/data/geodata/WDPA_WDOECM_Jun2022_Public_AF_shp_2.zip")
prot3 <- crop(readOGR("/work/idumville/pangolins/RAD/data/geodata/WDPA_WDOECM_Jun2022_Public_AF_shp-polygons.shp"), Cam)


prot1 <- sf::st_as_sf(prot1,  crs = 4326)
prot2 <- sf::st_as_sf(prot2,  crs = 4326)
prot3 <- sf::st_as_sf(prot3,  crs = 4326)

pnts_sf <- sf::st_as_sf(pred, coords = c("kd_x", "kd_y"), crs = 4326)
pnts_trans <- st_transform(pnts_sf, "EPSG:32632") 
prot1_trans <- st_transform(prot1, "EPSG:32632") 
prot2_trans <- st_transform(prot2, "EPSG:32632")  
prot3_trans <- st_transform(prot3, "EPSG:32632")   
res <- sf::st_join(pnts_trans, prot1_trans)
res <- sf::st_join(res, prot2_trans)
res <- sf::st_join(res, prot3_trans)
print(res)
resdf <- as.data.frame(res)
#pred$park <- resdf[c(which(resdf$NAME.x != "NA")), c("ORIG_NAME.x")] 

parks <- resdf[c(which(resdf$NAME.x != "NA")), c("ORIG_NAME.x")]
#
parks <- resdf[, c("sampleID", "ORIG_NAME.x")]

prot <- merge(as(prot1_trans, "Spatial"), as(prot2_trans, "Spatial"),  as(prot3_trans, "Spatial"))

protdist <- apply(gDistance(as(pnts_trans, "Spatial"), prot,byid=TRUE),2,min)

pred$parkdist <- protdist


##cameroon administration
 unzip("administrative-areas-latitude-and-longitude-coordinates-camadminshp.zip")

#camadm <- crop(read_OGR("/work/idumville/pangolins/RAD/data/geodata/CAM.shp"), Cam)
camadm <- read_sf("/work/idumville/pangolins/RAD/data/geodata/CAM.shp", crs = 4326)

pnts_sf <- sf::st_as_sf(pred, coords = c("kd_x", "kd_y"), crs = 4326)
pnts_trans <- st_transform(pnts_sf, "EPSG:32632") 
camadm_trans <- st_transform(camadm, "EPSG:32632") 
res <- sf::st_join(pnts_trans, camadm_trans)

resdf <- as.data.frame(res)
area <- resdf[,c("sampleID", "MMT_ID", "SHORT__FRM", "LONG_FRM", "ADM0", "ADM1", "ADM2")]
pred <- merge(pred, area, by="sampleID", all=T)


#get market distance for 90 and 00% RECALCULATING FOR INDV TOO BECUASE DIFFERENT

marketcoord<-read.delim("/work/idumville/pangolins/RAD/results/locator/concat.allWCA/allWCAknownmetadata.txt")
#pred <-  read.delim("../../results/locator/indvidualcentroids.txt")

temp <- merge(marketcoord, pred, by="sampleID")
temp <- temp[,c(1,3,2,6,7)] #swap x and y for market coordinates

pred$dist.km <- sapply(1:nrow(temp),function(i)
                spDistsN1(as.matrix(temp[i,2:3]),as.matrix(temp[i,4:5]),longlat=T))
                
################################
pred <- merge(pred, parks, by="sampleID", all=T) #whichparks are they in

#write.table(pred, "/work/idumville/pangolins/RAD/results/locator/indvidualpredictionswgeo.txt", quote=FALSE, row.names=FALSE, sep="\t")
#write.table(pred, "/work/idumville/pangolins/RAD/results/locator/00predictionswgeo.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(pred, "/work/idumville/pangolins/RAD/results/locator/90predictionswgeo.txt", quote=FALSE, row.names=FALSE, sep="\t")


#########################################
######## measuring distacne between reps ##
#######################################
#repeat for indv, 00 and 90

repcen=read.table("/work/idumville/pangolins/RAD/results/locator/indvrepdist.txt")

repcen$dist.km <- sapply(1:nrow(repcen),function(i)
                spDistsN1(as.matrix(repcen[i,2:3]),as.matrix(repcen[i,5:6]),longlat=T))
write.table(repcen, "/work/idumville/pangolins/RAD/results/locator/indvrepdist.txt", quote=FALSE, row.names=FALSE, sep="\t")



########## I can't get this to work with combn, had rows on the same line for same indv in old frame

combn(x, m, FUN = NULL, simplify = TRUE, …)

combs <- combn(seq_len(nrow(repcen)), 2)
Comb <- spDistsN1(as.matrix(repcen[combs[1,], ]), as.numeric(repcen[combs[2,], ]), longlat=T)
rownames(Comb) <- apply(combn(rownames(Dataset), 2), 2, paste, collapse = " ")


recpen <- as.matrix(repcen)
library(geosphere)
combn(seq_len(nrow(repcen)), 2, FUN= function(x) spDistsN1(repcen[x[1],], repcen[x[2],]))
#distGeo(repcen[x[2],], repcen[x[3],]))



