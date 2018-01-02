# Separate code to read and process LiDAR along Brown's Creek
# Load libraries
library(ggplot2)
library(data.table)
library(reshape2)
library(tidyr)
library(corrplot)
library(rgdal)
library(maptools)
library(lattice)
library(sp)
library(sf)  #Simple Features
library(raster)
library(RColorBrewer)
library(wesanderson)  #another color option
library(classInt)
library(lidR)

#Read in spatial data
#Read two LAZ tiles
#las1 <- readLAS("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/las/3542-30-37.laz")
#las2 <- readLAS("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/las/3542-30-38.laz")
#read in spatial data from shapefiles 
shp.boundary <- st_read("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/Boundary2.shp")  #Study Area Bounndary
proj4 <- "+proj=utm +zone=15"

#Write las files to folder
#writeLAS(las1,"H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/las_R_output/3542-30-37.las")
#writeLAS(las2,"H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/las_R_output/3542-30-38.las")

#Plot to check that they look good
#plot(lidar1)
#plot(lidar2)

#Build catalogue of las tiles
ctg <- catalog("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/test/")
newctg = catalog_reshape(ctg, 6000, "H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/catalog", "lidar_2")
#plot(ctg, proj4)

#create terrain (bare earth) model (Digital Elevation Model)
#dem1 <- grid_terrain(las1, method = "knnidw", k = 6)

#Create a canopy surface model with 4m^2 cells
grid_canopy(ctg) %>% plot
#Mean height with 400 m^2 cells
#grid_metrics(ctg, mean(Z)) %>% plot
#With multiple metrics
#metrics = grid_metrics(ctg, .stdmetrics_z)
plot(metrics)

#merge and clip the tiles into 1 las file



#determine max height of veg within 5 m buffer of creek




#determine density of vegetation canopy within 5 m buffer of creek




#outputs to use in other code



# Outliers of intensity breaks the color range. Use the trim parameter.
#plot(lidar, color = "Intensity", colorPalette = heat.colors(50))
#plot(lidar, color = "Intensity", colorPalette = heat.colors(50), trim = 0.99)