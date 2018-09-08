#Brown's Creek Riparian Shading Study September 6, 2018
#by Olivia Sparrow
#Civil Engineering Master's Student, University of Minnesota-Twin Cities
#Water Resources Engineer, Emmons & Olivier Resources

#Description of code: This code reads in the shade estimated using two methods:
#    (1) LiDAR (collected in 2011, analyzed by Herb & Correll in 2016) processed in Arcmap radiation tool
#    (2) hemispherical photos (collected by me in 2017) processed in WinSCANOPY
#
#The code also reads in physical characteristics of the monitoring sites. The code manipulates the data
#to assess potential predictive models of shade.


# Section 1: Notes ---------------------------------------------------


#left and right banks identified looking downstream


# Section 2: Libraries ----------------------------------------------------


#load libraries
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
library(classInt)
library(lidR)
library(plyr)
library(wesanderson)
library(lemon)

# Section 3: Read and format data -----------------------------------------

#Read in and format the field data
df <- read.csv('TidyTransectData.csv', header = TRUE)  #Read full table of data collected
colnames(df)[colnames(df) =='reach.id'] <- 'Reach'
df <- df[!(df$Reach==7),]   #remove reach 7 from analysis since this reach was restored in 2012 and cannot be used to
                            #relate LiDAR (2011) and hemi photo (2017) results 
df$date <- as.character(df$date)  #convert date to character string
df$time.m <- as.character(df$time.m)  #convert time of middle photo to character string
df$time.l <- as.character(df$time.l)  #convert time of left photo to character string
df$time.r <- as.character(df$time.r)  #convert time of right photo to character string
df$time.m <- as.POSIXct(paste(df$date, df$time.m), format = '%m/%d/%Y %H:%M')  #add date of middle photo to timestamp
df$time.l <- as.POSIXct(paste(df$date, df$time.l), format = '%m/%d/%Y %H:%M')  #add date of left photo to timestamp
df$time.r <- as.POSIXct(paste(df$date, df$time.r), format = '%m/%d/%Y %H:%M')  #add date of right photo to timestamp
df$date <- as.Date(df$date, format = '%m/%d/%Y')  #convert date column from factor to date
#add the necessary metadata to plot the lat/lon coordinates
df$geodetic.da <- "NAD83"
df$utm.zone <- "15N"

#read in spatial data from shapefiles 
shp.boundary <- st_read("D:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/Boundary2.shp")  #Study Area Bounndary
shp.creek <- st_read("D:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/Creek.shp")  #Creek line
shp.shade <- st_read("D:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/Raw LiDAR shade results/Model_SegsBuffer2.shp")  #Thermal Study Shade analysis
shp.shade.dem <- st_read("D:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/out_global_radiation2.shp")  #Topographic shade
#Height of vegetation above river at transects from DSM

#add northings and eastings to df since coordinates are only as lat and long format
utm15nCRS <- st_crs(shp.shade)  #save geospatial metadata
cord.dec <- SpatialPoints(cbind(df$long.dec.deg, df$lat.dec.deg), 
                          proj4string = CRS("+proj=longlat"))  #save existing coodinates in lat-long system
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:26915"))
df.cord.utm <- as.data.frame(cord.UTM)  #dataframe of coordinates of all points monitored in study
df$easting <- df.cord.utm$coords.x1
df$northing <- df.cord.utm$coords.x2
#write.csv(df, "data_table_with_coords") #write df to csv for use in ArcGIS

# Section 4: Spatial Data Analysis --------------------------------------------

#TODO add back in the analysis of the calibrated Lidar shade from the thermal study

#plot the project boundary and creek location (Note that coordinates are in Northing and Easting)
plot(st_geometry(shp.boundary), border = "black", lwd = 2, main = "Study Area Boundary")  #plots boundary of study area
plot(st_geometry(shp.creek), add = TRUE, col = "blue")

#calculate total growing season shade from 2011 LiDAR analysis
months <- c("MayMN", "JunMN", "JulMN", "AugMN", "SepMN")  #define headings of columns used in calculations below
df.shp.shade <- as.data.frame(shp.shade)  #create a dataframe version of the shade shapefile attributes
df.shp.shade$total.rad <- rowSums(df.shp.shade[, c(months)])  #sum of total radiation reaching each transect over the growing season
df.shp.shade$max.rad <- sum(sapply(df.shp.shade[, c(months)], max))  #sum of maximum radiation each month (i.e. "above canopy" radiation)
df.shp.shade$Shade <- (1 - df.shp.shade$total.rad / df.shp.shade$max.rad) * 100  #calculate total shade at each transect
shp.shade$Shade <- df.shp.shade$Shade  #save shade back to shapefile
shp.shade <- shp.shade[shp.boundary,]  #clip to study area boundary

#plot the shade estimated by GIS analysis of 2011 LiDAR (Note that coordinates are in Northing and Easting)
pal2 <- rev(wes_palette(n = 5, name = "Zissou1"))
br <- c(0, 20, 40, 60, 80, 100)
offs <- 0.0000001 
br[1] <- br[1] - offs 
br[length(br)] <- br[length(br)] + offs 
shade_cuts <- cut(shp.shade$Shade, br)
#plot shade estimated using ArcGIS Solar Radiation tool and LiDAR
plot(shp.shade["Shade"], lwd = 0.2, axes = TRUE, col = pal2[as.numeric(shade_cuts)], border = pal2[as.numeric(shade_cuts)],
     xlab = "Easting", ylab = "Northing", main = "Shade from GIS Analysis of 2011 LiDAR", cex.axis = 0.6)  
legend("topright", legend = paste("<", round(br[-1])), col = pal2, lty = 1, lwd = 2, cex = .3)

#add transect locations to the above plot (Note that coordinates are in lat and long)
pts.transects <- st_as_sf(df, coords = c("easting", "northing"), crs = utm15nCRS)  #convert field data to sf for map
plot(shp.shade["Shade"], lwd = 2, axes = TRUE, col = pal2[as.numeric(shade_cuts)],
     border = pal2[as.numeric(shade_cuts)], xlab = "Easting", ylab = "Northing", cex.axis = 0.6,
     main = "Shade from GIS Analysis of 2011 LiDAR & Riparian Shading Study Transects")  #plot shade estimated using ArcGIS Solar Radiation tool and LiDAR
legend("topright", legend = c("Transects", paste("<", round(br[-1]))), pch = c(1, NA, NA, NA, NA, NA, NA), 
       lty = c(NA, 1, 1, 1, 1, 1, 1), col = c("black", pal2), lwd = 2, cex = .3)
plot(st_geometry(pts.transects), add = TRUE)  # add transect locations to plot

#determine and save the shade estimated using lidar into the df
dist <- st_distance(pts.transects, shp.shade)  #calculate distance between each creek segment and transect point
i <- as.matrix(apply(dist, 1, which.min))  #find the row index in dist with the shortest distance to each point
ss <- shp.shade$Shade[i]
pts.transects$shade.lidar.leafoff.dsm <- ss
df[, "shade.lidar.leafoff.dsm"] <- ss
shade_cuts <- cut(pts.transects$shade.lidar.leafoff.dsm, br)
plot(pts.transects["shade.lidar.leafoff.dsm"], axes = TRUE, col = pal2[as.numeric(shade_cuts)],
     xlab = "Easting", ylab = "Northing", cex.axis = 0.6,
     main = "Shade from GIS Analysis of 2011 LiDAR at Riparian Shading Study Transects")
legend("topright", legend = paste("<", round(br[-1])), col = pal2, pch = 1, lwd = 2, cex = .3)
df$shade.lidar.leafoff.dsm <- df$shade.lidar.leafoff.dsm / 100  #convert to fraction

#TODO: the above uses the lidar results AFTER they were calibrated by Bill Herb - check his email to see if 
#I need to use original results. They were also averaged over 20+ foot segments of the stream, but Bill did not recommend
#going to original results because they were so noisy. As such, it might be more appropriate to compare the average shade
#across each stream segment resulting from the lidar and hemiphoto riparian shade analysis instead of at specific points.

#Read the vegetation height above river (HAR) raster and populate the dataframe with this new information
shp.pt.circ <- st_buffer(pts.transects, dist = 10)  #create 27m buffer circles around each data point
har.left <- raster("D:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/rasters/Canopy_HAR_Left1.tif")  #read left bank HAR
har.right <- raster("D:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/rasters/Canopy_HAR_Right1.tif")  #read left bank HAR
buf.elev.l <- raster::extract(har.left, shp.pt.circ, fun = max, sp = TRUE, na.rm = TRUE)  #extract max HAR into buffers around data pts
buf.elev.r <- raster::extract(har.right, shp.pt.circ, fun = max, sp = TRUE, na.rm = TRUE)  #extract max HAR into buffers around data pts
df$HAR.l.m <- buf.elev.l$Canopy_HAR_Left1
df$HAR.r.m <- buf.elev.r$Canopy_HAR_Right1
df$veg.height.max.l.m <- apply(df[,c("HAR.l.m","herb.height.max.l.m")], 1, #save maximum veg height from field
                             function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))  #and LiDAR data
df$veg.height.max.r.m <- apply(df[,c("HAR.r.m","herb.height.max.r.m")], 1, #save maximum veg height from field
                               function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))  #and LiDAR data
  #max function would return NA values without the code at the end of the above command

# Section 5: Subset Dataframes -------------------------------------------

#subset the data for use in analysis of stage-shade curves (other subsetting done after correction to common height)
df.stage <- subset(df, purpose == 'Stage')  #Create a smaller table of just the staging observations

#Break the stage and height dataframes into separate objects so they can be melted and then merged later
df.stage.height <- data.frame("Reach" = df.stage$Reach, "Transect" = df.stage$transect.no,
                              "Left" = df.stage$lens.height.l, "Middle" = df.stage$lens.height.m,  
                              "Right" = df.stage$lens.height.r, "Azimuth" = df.stage$gen.strm.azimuth,
                              "Vegetation" = df.stage$veg.type.l)
df.stage.shade <- data.frame("Reach" = df.stage$Reach, "Transect" = df.stage$transect.no,
                              "Left" = df.stage$shade.l, "Middle" = df.stage$shade.m,  
                              "Right" = df.stage$shade.r, "Azimuth" = df.stage$gen.strm.azimuth,
                             "Vegetation" = df.stage$veg.type.l)
#melt and rename stage (height) and shade data frames by position
df.stage.height.m <- melt(df.stage.height, id.vars = c("Reach", "Transect", "Azimuth", "Vegetation"), 
                          measure.vars = c("Left", "Middle", "Right"))
colnames(df.stage.height.m) <- c("Reach", "Transect", "Azimuth", "Vegetation", "Position", "lens.height")
df.stage.shade.m <- melt(df.stage.shade, id.vars = c("Reach", "Transect", "Azimuth", "Vegetation"), 
                          measure.vars = c("Left", "Middle", "Right"))
colnames(df.stage.shade.m) <- c("Reach", "Transect", "Azimuth", "Vegetation", "Position", "shade")
#merge staging dataframes
df.stage.m <- data.frame(df.stage.height.m, "shade" = df.stage.shade.m$shade)
setorder(df.stage.m, lens.height)
df.stage.m <- mutate(df.stage.m, Vegetation = revalue(Vegetation, c("F" = "Woody", "G" = "Grassy", "M" = "Mixed")))
#subset staging df by transect and position
df.stage.m.36L <- subset(df.stage.m, Transect == 6 & Position == "Left")
df.stage.m.36M <- subset(df.stage.m, Transect == 6 & Position == "Middle")
df.stage.m.36R <- subset(df.stage.m, Transect == 6 & Position == "Right")
df.stage.m.39L <- subset(df.stage.m, Transect == 9 & Position == "Left")
df.stage.m.39M <- subset(df.stage.m, Transect == 9 & Position == "Middle")
df.stage.m.39R <- subset(df.stage.m, Transect == 9 & Position == "Right")
#write each to curve CSV
write.csv(df.stage.m, "data_table_staging_tables")


# Section 6: Regression of Stage-Shade Curves and Correction --------------

#Remove outlier point on Reach 3, Transect 9, Left at highest elevation above stream
df.stage.m <- df.stage.m[!(df.stage.m$Transect == 9 & df.stage.m$lens.height > 0.7 & df.stage.m$Position == "Left"),]

#########################
######REPORT FIGURE######
#Plot shade-stage curves for 4 reference transects
ggplot(df.stage.m, aes(lens.height, shade * 100, colour = Position)) + geom_point() + geom_line() +
  labs(x = "Height of Lens Above Water (m)", y = "HP-Based Shade, Growing Season Average (%)") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  facet_wrap(~ Vegetation + Azimuth,  labeller = label_both) +
  theme(axis.text = element_text(size = 18)) + 
  theme_grey(base_size = 20) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) + theme(legend.position = "bottom")

#Define standard x value (xs) and interpolate corresponding ys value
#determine standard reference values for x and y
xs <- 0.2
rr <- c(3, 3, 3, 3, 3, 3, 2, 2, 2, 4, 4, 4)
tt <- c(6, 6, 6, 9, 9, 9, 0, 0, 0, 5, 5, 5)
pp <- c("Left", "Middle", "Right", "Left", "Middle", "Right", "Left", "Middle", "Right", "Left", "Middle", "Right")
x1 <- as.numeric(1)
x2 <- as.numeric(1)
y1 <- as.numeric(1)
y2 <- as.numeric(1)
ys <- as.numeric(1)
df.ys <- data.frame(rr, tt, pp, x1, x2, y1, y2, xs, ys)
for(i in 1:length(rr))  #linear interpolation variables (x1, x2, y1, and y2) 
{
  r <- df.ys[i, 1]  #current reach
  t <- df.ys[i, 2]  #current transect
  p <- df.ys[i, 3]  #current position (left, middle, right)
  df.ys[i, 4] <- min(df.stage.m$lens.height[df.stage.m$Reach == r & df.stage.m$Transect == t &
                                     df.stage.m$Position == p])
  df.ys[i, 5] <- min(df.stage.m$lens.height[df.stage.m$Reach == r & df.stage.m$Transect == t &
                                     df.stage.m$Position == p & df.stage.m$lens.height > df.ys[i, 4]])
  df.ys[i, 6] <- df.stage.m$shade[df.stage.m$Reach == r & df.stage.m$Transect == t &
                                     df.stage.m$Position == p & df.stage.m$lens.height == df.ys[i, 4]]
  df.ys[i, 7] <- df.stage.m$shade[df.stage.m$Reach == r & df.stage.m$Transect == t &
                           df.stage.m$Position == p & df.stage.m$lens.height == df.ys[i, 5]]

}
#interpolate reference value for shade from each stage-shade curve
df.ys$ys <- (df.ys$xs - df.ys$x1) * (df.ys$y2 - df.ys$y1) / (df.ys$x2 - df.ys$x1) + df.ys$y1  
names(df.ys)[1] <- "Reach"
names(df.ys)[2] <- "Transect"
names(df.ys)[3] <- "Position"
#save reference values to staging df
df.stage.mer <- merge(df.stage.m, df.ys, by = c("Reach", "Transect","Position"), all.x = TRUE)
#normalize the height and shade using reference values
df.stage.mer$xstar1 <- df.stage.mer$lens.height / df.stage.mer$xs  #standardize lens height above stream
df.stage.mer$ystar1 <- df.stage.mer$shade / df.stage.mer$ys  #standardize shade
#plot normalized shade-stage curves for all reference transects
ggplot(df.stage.mer, aes(xstar1, ystar1, colour = Position)) + geom_point() + geom_line() +
  labs(title = "1b: Standardized Stage-Shade Curves (WinSCANOPY)", x = "Normalized Lens Height Above Water",
       y = "Normalized Average Shade Over Growing Season") +
  facet_wrap(~ Reach + Transect,  labeller = label_both)
df.stage.mer$line <- paste(df.stage.mer$Reach, df.stage.mer$Transect, df.stage.mer$Position)

#plot normalized shade-stage curves in one plot
ggplot(df.stage.mer, aes(xstar1, ystar1, colour = Position, group = line)) + geom_point(aes(size = 1, shape = Vegetation)) + geom_line() +  #without wrapping to see consistency
  labs(x = "Normalized Lens Height Above Water", y = "Normalized HP-Based Shade") +
  theme(axis.text = element_text(size = 18)) + 
  theme_grey(base_size = 20) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) + theme(legend.position = "bottom")

#Developing regression models to use in correcting shade in the transects df
#divide the df into three based on shade since they look like they have three separate relationships
df.stage.mer1 <- subset(df.stage.mer, shade < 0.5 & Position != "Left")  #df for developing model for adjusting 
                                                                        #low-shaded transects in middle and right positions
df.stage.mer2 <- subset(df.stage.mer, shade >= 0.5)  #df for high-shaded areas
df.stage.mer3 <- subset(df.stage.mer, Position == "Left" & Transect == 6)  #df for the unique left position curve
#note: df.stage.mer2 and df.stage.mer3 are not used beyond this point because forested areas did not see a lot of change vs. height
    #for the left bank, there may be significant change from heightened elevations, however these were not common in grassy reaches.
    #The left position in transect 3-6 was an outlier. 
    #Hereafter, adjustments are only made for middle and right bank shading in low-shaded reaches.

#define independent variables for plotting and modeling
hh1 <- df.stage.mer1$xstar1
ss1 <- df.stage.mer1$ystar1

# Nonlinear Regression: shade as function of lens height above stream
# Solution using Downhill Simplex method
data <- data.frame(cbind(c(1:20), hh1, ss1))
names(data) <- c("obs", "norm.lens.height", "norm.shade")  # Rename column headings
data  # Print original data in file
summary(data)  # Basic Statistics
x = data[, 2]  # Plot using second and third columns for lens height and shade
y = data[, 3]
plot(x, y, xlab = 'Normalized height above stream', ylab = 'Normalized shade') 
fexample1 <- function(b_parm) {  # Define minimization function to use in solving for the parameters
  x = data[, 2]
  y = data[, 3]
  bo = b_parm[1]
  b1 = b_parm[2]
  sum( ( y - bo * exp(b1 * x) ) ^ 2 )
}                     
# Determine optimal values using downhill simplex method of Nelder and Mead
# Use initial values of bo = 10, and b1 = -1
nonlinear = optim(c(10, -1), fexample1)  # find the optimal parameters for the regression line
nonlinear$par  # Best fit parameter values
nonlinear$value  # Residual sum of squares
fexample2 <- function(b_parm, x) {  # Define exponential regression function
  bo = b_parm[1]
  b1 = b_parm[2]
  bo * exp(b1 * x)
}
fexample2(nonlinear$par, 2)  # Illustration of computation of predicted value
# Compute ANOVA values (anova(nonlinear) not available for nls)
SSTO = sum((data$Norm.shade - mean(data$norm.shade)) ^ 2)
SSE = nonlinear$value
c_d <- 1 - SSE/SSTO
SSR = sum((fexample2(nonlinear$par, x)
           - mean(data$norm.shade)) ^ 2)
sprintf('SSTO=%.1f, SSE=%.1f, SSR=%.1f, c_d=%.3f', SSTO, SSE, SSR, c_d)
# Graph solution		
par(cex = 2)
plot(norm.shade ~ norm.lens.height, data, pch = 2,
     xlab = 'Normalized lens height above stream', ylab = 'Normalized shade')
r = range(data$norm.lens.height)
d = seq(r[1], r[2], length = 100)
lines(d, fexample2(nonlinear$par, d))
legend(2.8, 1.1, c('obs', 'fit'), pch = c(2, NA), lty = c(NA, 1))

# Create function to correct measured x and y to a common height (xc) above stream using the above exponential regression 
fcorr <- function(b_parm, xm, ym, xs = 0.2, xc = 0.1){
  bo = b_parm[1]
  b1 = b_parm[2]
  ym * exp(b1 * (xc - xm) / xs)
}

# Update data frames to save the old Winscanopy shade (at varying lens heights)
df$shade.varh.m <- df$shade.m
df$shade.varh.l <- df$shade.l
df$shade.varh.r <- df$shade.r
df$shade.varh.avg <- df$shade.avg

# Apply the function to correct height-varying shade
# Correction was only applied to middle and right positions where original shade was less than 50%
for(i in 1:nrow(df)){  #linear interpolation variables (x1, x2, y1, and y2) 
  if (df[i, which(colnames(df) == "shade.varh.m")] < 0.5) {
    df[i, which(colnames(df) == "shade.m")] <- fcorr(b_parm = nonlinear$par, 
                                                   xm = df[i, which(colnames(df) == "lens.height.m")], 
                                                   ym = df[i, which(colnames(df) == "shade.varh.m")])
  }
  
  if (df[i, which(colnames(df) == "shade.varh.r")] < 0.5) {
    df[i, which(colnames(df) == "shade.r")] <- fcorr(b_parm = nonlinear$par, 
                                                     xm = df[i, which(colnames(df) == "lens.height.r")], 
                                                     ym = df[i, which(colnames(df) == "shade.varh.r")])
  }
  
  df[i, which(colnames(df) == "shade.avg")] <- (df[i, which(colnames(df) == "shade.m")] + df[i, which(colnames(df) == "shade.l")] +
                                                  df[i, which(colnames(df) == "shade.r")]) / 3
}
#Subset dataframe for use in the subsequent sections
df.transect <- subset(df, purpose == 'Transect')  #Create a smaller table of just the transect observations
df.transect.num <- df.transect[, sapply(df.transect, is.numeric)]  #Create df of only numeric values for correlation analysis
df.transect.grass <- subset(df.transect, veg.code.avg == 0)  #Numerical observations at grassy transects
df.test <- subset(df,purpose == 'Test')  #Create a smaller table of just the two sites to be used to test the regression

#Create long dataframes for plotting results for each transect and stage-shade curves
df.transect.m <- melt(df.transect, id.vars = c("Reach", "transect.no"), 
                      measure.vars = c("shade.l", "shade.m", "shade.r", "shade.avg"))
colnames(df.transect.m)[colnames(df.transect.m) =='variable'] <- 'Position'  #change column heading for position across transect
#rename positions for graphing later
df.transect.m <- mutate(df.transect.m, Position = revalue(Position, c("shade.l" = "Left", "shade.m" = "Middle", "shade.r" = "Right", "shade.avg" = "Average")))
df.transect.varh.m <- melt(df.transect, id.vars = c("Reach", "transect.no"), 
                      measure.vars = c("shade.varh.l", "shade.varh.m", "shade.varh.r", "shade.varh.avg"))
#export shade results for use in stream temperature model shade input
vars <- c("date", "Reach", "transect.no", "station.m", "shade.lidar.leafoff.dsm", "shade.avg")  #define variables to keep in exported df
dfclip <- df[vars]  #create clipped version of df
write.csv(dfclip, file = "data_table_winscanopy_transect_shade")
write.csv(df, "data_table_all.csv")  #Write a copy of the full data table to a csv file


# Section 7: Correlation Matrix -------------------------------------------


#assess correlation between independent & dependent variables
df.transect.num$total.site.factor.m <- NULL  #delete total site factor columns
df.transect.num$total.site.factor.l <- NULL  #delete total site factor columns
df.transect.num$total.site.factor.r <- NULL  #delete total site factor columns
df.transect.num$total.site.factor.avg <- NULL  #delete total site factor columns
df.transect.corr <- cor(df.transect.num, use = "na.or.complete")
corrplot.mixed(df.transect.corr, upper = "color", number.cex = .4, number.digits = 1, lower.col = "black", 
               tl.cex = 0.5, tl.pos = "lt", tl.col = "black", cl.cex = 0.5, 
               title = "2a: Correlation coefficients of numeric variables (All Transects)", mar=c(0,0,1,0))
df.transect.num.grass <- subset(df.transect.num, veg.code.avg == 0)  #Numerical observations at grassy transects
df.transect.num.grass$veg.code.avg <- NULL
df.transect.num.grass$veg.code.l <- NULL
df.transect.num.grass$veg.code.r <- NULL
df.transect.grass.corr <- cor(df.transect.num.grass, use = "na.or.complete")
corrplot.mixed(df.transect.grass.corr, upper = "color", number.cex = .4, number.digits = 1, lower.col = "black", 
               tl.cex = 0.5, tl.pos = "lt", tl.col = "black", cl.cex = 0.5, 
               title = "2b: Correlation coefficients of numeric variables (Grassy transects)", mar=c(0,0,1,0))


# Section 8: Plot Observed Data and Summary Statistics --------------------

#how many transects are entirely grassy
length(which(df.transect$veg.type.l == "G", df.transect$veg.type.r == "G"))

#########################
######REPORT FIGURE######
#Plot shade for each reach and transect and position (after lens height correction)
ggplot(df.transect.m, aes(transect.no, 100 * value, colour = Position)) + geom_point() + geom_line() +
  labs(x = "Transect", y = "HP-Based Shade, Growing Season Avg (%)", caption = "after correction for lens height above stream") +
  scale_x_continuous(breaks = seq(0, 11, 1)) + scale_y_continuous(breaks = seq(0, 100, 20)) +
  facet_wrap(~ Reach,  labeller = label_both, ncol = 2) + theme(axis.text = element_text(size = 18)) + 
  theme_grey(base_size = 20) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) + theme(legend.position = "bottom")

summary(df.transect)  #summary statistics of the transects (not stage or test locations) used in reported statistics

#Plot shade for each reach and transect and position (before lens height correction)
# ggplot(df.transect.varh.m, aes(transect.no, 100 * value, colour = variable)) + geom_point() + geom_line() +
#   labs(title = "1b. Estimated Shade from WinSCANOPY Simulation", x = "Transect #", 
#        y = "Average % Shade Over Growing Season", caption = "before correction for lens height above stream") +
#   scale_x_continuous(breaks = seq(0, 11, 1)) + scale_y_continuous(breaks = seq(0, 100, 20)) +
#   facet_wrap(~ Reach,  labeller = label_both)

#Plot average shade at each transect and reach
ggplot(df.transect, aes(transect.no, shade.avg * 100, group = factor(Reach), colour = factor(Reach))) +
  geom_point() + geom_line() + scale_x_continuous(breaks = seq(0, 11, 1)) + scale_y_continuous(breaks = seq(0, 100, 20)) +
  labs(title = "2a. Average Shade at Each Transect from WinSCANOPY Simulation",
       x = "Transect #", y = "Average % Shade Over Growing Season") +
  facet_wrap(~ Reach,  labeller = label_both)

#########################
######REPORT FIGURE######
#Plot average shade at each transect (grouped by reach) - comparing hemiphotos to LiDAR
df.transect.shade.m <- melt(df.transect, id.vars = c("Reach", "transect.no"), 
                      measure.vars = c("shade.avg", "shade.lidar.leafoff.dsm"))
colnames(df.transect.shade.m)[3] <- "Method"  #update to reflect that variable column is the method of estimating shade
levels(df.transect.shade.m$Method) <- c(levels(df.transect.shade.m$Method), "HP", "LiDAR")
df.transect.shade.m$Method[df.transect.shade.m$Method == "shade.avg"] <- "HP"
df.transect.shade.m$Method[df.transect.shade.m$Method == "shade.lidar.leafoff.dsm"] <- "LiDAR"
ggplot(df.transect.shade.m, aes(transect.no, 100 * value, group = factor(Reach), colour = Method)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 11, 1)) + scale_y_continuous(breaks = seq(0, 100, 20)) +
  labs(x = "Transect", y = "Shade, Growing Season Avg (%)") +
  facet_wrap(~ Reach,  labeller = label_both) + theme(axis.text = element_text(size = 18)) + 
  theme_grey(base_size = 20) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) + theme(legend.position = "bottom")

#Plot lidar-derived shade vs. Winscanopy derived shade *by transect*
ggplot(df.transect, aes(x = shade.lidar.leafoff.dsm, y = shade.avg * 100, color = veg.code.avg)) + geom_point() +
  scale_colour_gradientn(colours = rainbow(4))
ggplot(df.transect, aes(x = veg.code.avg , y = openness.avg)) + geom_point() +
  scale_colour_gradientn(colours = rainbow(4))

#calculate and plot average LiDAR v WinsCANOPY shade *by reach*
df.transect.agg <- aggregate(df.transect[, c('shade.avg', 'shade.lidar.leafoff.dsm')], list(df.transect[, c('Reach')]), mean)
ggplot(df.transect.agg, aes(x = shade.lidar.leafoff.dsm * 100, y = shade.avg * 100)) + geom_point()
xx <- df.transect.agg$shade.lidar.leafoff.dsm
yy <- df.transect.agg$shade.avg
write.csv(df.transect.agg, "data_table_avg_reach_shade") 

#Linear regression model (Lidar vs Winscanopy shade)
mod.reach <- lm(shade.avg ~ shade.lidar.leafoff.dsm, df.transect.agg)
anova(mod.reach)  #Print ANOVA results
summary(mod.reach)  #Print coefficients
df.reach <- data.frame(x = xx, y = predict(mod.reach, data.frame(shade.lidar.leafoff.dsm = xx)), 
                   sres = rstandard(mod.reach))  #df of standard residuals and predicted values
ggplot(df.reach, aes(x = x, y = sres)) + geom_point() +  #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Shade from Lidar") + ylab("Standardized Residuals") +
  ggtitle("Model: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect.agg, aes(x = shade.lidar.leafoff.dsm, y = shade.avg)) +
  geom_line(data = df.reach, aes(x = x, y = y)) + ylab("Shade from Winscanopy") +
  xlab("Shade from Lidar") +
  ggtitle("Model: Observed & Predicted Values")

#Plot observed values vs. average transect shade (grassy only)
df.transect.num.grass.m <- melt(df.transect.num.grass, id = c('shade.avg'))
ggplot(data = df.transect.num.grass.m, aes(x = value, y = shade.avg * 100, colour = variable)) +
  geom_point() +  xlab("Independent Variables") + ylab("Average Transect Shade Over Growing Season (%)") +
  ggtitle("6a. Independent Variables vs. Avg Transect Shade (Grassy Transects)") + 
  facet_wrap( ~ variable, scales="free_x") + theme(legend.position = "none")

#Plot observed values vs. average transect shade (all vegetation types)
df.transect.num.m <- melt(df.transect.num, id.vars = c('shade.avg'))
ggplot(data = df.transect.num.m, aes(x = value, y = shade.avg * 100, colour = variable)) +
  geom_point() +  xlab("Independent Variables") + ylab("Average Transect Shade Over Growing Season (%)") +
  ggtitle("6b. Independent Variables vs. Avg Transect Shade (All Transects)") + scale_y_continuous(breaks = seq(0, 100, 20)) +
  facet_wrap( ~ variable, scales="free_x") + theme(legend.position = "none")

#########################
######REPORT FIGURE######
#Plot shade based on vegetation type
df.transect.m.shadeveg <- melt(df.transect, id.vars = c("Reach", "transect.no", "veg.code.avg"), 
                                   measure.vars = c("shade.avg"))
df.transect.m.shadeveg$veg.code.avg <- factor(df.transect.m.shadeveg$veg.code.avg * 100)
ggplot(df.transect.m.shadeveg, aes(veg.code.avg, value * 100)) + geom_boxplot(outlier.size = 3) + 
  scale_y_continuous(breaks = seq(0, 100, 20)) + xlab("Riparian Vegetation Compisition: % Woody") +
  theme_grey(base_size = 20) + ylab("Average Transect Shade Over Growing Season (%)") + 
  theme(axis.text = element_text(size = 18)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

#Plot vegetation height based on vegetation type
df.transect.num.m.height <- subset(df.transect.num.m, variable == "herb.height.max.l.m" | 
                                           variable == "herb.height.max.r.m" | 
                                           variable == "herb.height.mode.l.m" | 
                                           variable == "herb.height.mode.r.m")
ggplot(df.transect.num.m.height, aes(variable, value)) + geom_boxplot() +
  ylab("Veg Height (m)") + ggtitle("7a. Veg Height (All Transects)")
#grassy only
df.transect.num.grass.m.height <- subset(df.transect.num.grass.m, variable == "herb.height.max.l.m" | 
                                           variable == "herb.height.max.r.m" | 
                                           variable == "herb.height.mode.l.m" | 
                                           variable == "herb.height.mode.r.m")
ggplot(df.transect.num.grass.m.height, aes(variable, value)) + geom_boxplot() +
  ylab("Veg Height (m)") + ggtitle("7b. Veg Height (Grassy Transects)")

#Plot wetted width based on vegetation type
df.transect.width.m <- melt(df.transect, id.vars = c("Reach", "transect.no", "veg.code.avg"), 
                              measure.vars = c("wetted.width.m"))
df.transect.width.m$veg.code.avg <- factor(df.transect.width.m$veg.code.avg)
ggplot(df.transect.width.m, aes(veg.code.avg, value)) + geom_boxplot() + xlab("Riparian Vegetation (Grass = 0...Forest = 1)") +
  ylab("Wetted Width (m)") + ggtitle("7c. Wetted Width (All Transects)")

#Plot stream azimuth based on vegetation type
df.transect.genazimuth.m <- melt(df.transect, id.vars = c("Reach", "transect.no", "veg.code.avg"), 
                            measure.vars = c("gen.strm.azimuth"))
df.transect.genazimuth.m$veg.code.avg <- factor(df.transect.genazimuth.m$veg.code.avg)
ggplot(df.transect.genazimuth.m, aes(veg.code.avg, value)) + geom_boxplot() + xlab("Riparian Vegetation (Grass = 0...Forest = 1)") +
  ylab("General Stream Azimuth (degrees from South)") + ggtitle("7d. General Stream Azimuth (All Transects)")

#Plot lens height based on vegetation type
df.transect.lens.m <- melt(df.transect, id.vars = c('Reach', 'transect.no', 'veg.code.avg'), 
                           measure.vars = c('lens.height.m', 'lens.height.l', 'lens.height.r'))
df.transect.lens.m$veg.code.avg <- factor(df.transect.lens.m$veg.code.avg)
means <- aggregate(lens.height ~ veg.code.avg, df.transect.lens.m, mean) #mean lens height for each veg type
ggplot(df.transect.lens.m, aes(veg.code.avg, value)) + geom_boxplot() + xlab("Riparian Vegetation (Grass = 0...Forest = 1)") +
  ylab("Lens Height Above Water (m)") + ggtitle("8. Height of Camera Lens Above Stream (All Transects)")

#Plot lens height based on shade at each position
colnames(df.transect.m) <- c("Reach", "transect.no", "position", "shade")
colnames(df.transect.lens.m) <- c("Reach", "transect.no", "veg.code.avg", "position", "lens.height")
levels(df.transect.m$position) <- c("left", "middle", "right", "average")
levels(df.transect.lens.m$position) <- c("middle", "left", "right")
df.transect.lens.shade.m <- merge(df.transect.lens.m, df.transect.m, 
                                  by = c("Reach", "transect.no","position"), all.x = TRUE)
ggplot(df.transect.lens.shade.m, aes(shade, lens.height, colour = veg.code.avg)) + geom_point() + xlab("Shade at each position (left, middle, right)") +
  ylab("Lens Height Above Water (m)") + ggtitle("9a. Height of Camera Lens Above Stream (All Transects, All Positiions)")
df.transect.lens.shade.m.left <- subset(df.transect.lens.shade.m, position == "left")
ggplot(df.transect.lens.shade.m.left, aes(shade, lens.height, colour = veg.code.avg)) + geom_point() + xlab("Shade at each position (left, middle, right)") +
  ylab("Lens Height Above Water (m)") + ggtitle("9b. Height of Camera Lens Above Stream (All Transects, Left Positions)")

# # Section 9: Regression Models --------------------------------------------
# 
# #Create vectors for the independent variables
# vv <- df.transect$veg.code.avg
# ww <- df.transect$wetted.width.m 
# tw <- df.transect$thalweg.depth.m
# wwtw <- df.transect$width.to.depth
# lhm <- df.transect$lens.height.m
# lhl <- df.transect$lens.height.l  
# lhr <- df.transect$lens.height.r
# wel <- df.transect$wetted.edge.l.m
# wer <- df.transect$wetted.edge.r.m
# hmaxavg <- df.transect$herb.height.max.avg.m 
# vhmaxr <- df.transect$veg.height.max.r.m
# azm <- df.transect$strm.azimuth
# azgen <- df.transect$gen.strm.azimuth
# vdens <- df.transect$veg.dens.avg
# hmodavg <- df.transect$herb.height.mode.avg.m
# 
# #Linear regression model (Model #1a: Veg Type vs. Shade in Middle)
# mod1a <- lm(shade.m ~ veg.code.avg, df.transect)
# anova(mod1a)  #Print ANOVA results
# summary(mod1a)  #Print coefficients
# df1a <- data.frame(x = vv, y = predict(mod1a, data.frame(veg.code.avg = vv)), 
#                    sres = rstandard(mod1a))  #df of standard residuals and predicted values
# ggplot(df1a, aes(x = x, y = sres)) + geom_point() +  #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Code for vegetation on both banks (0=grass; 1=forest)") + ylab("Standardized Residuals") +
#   ggtitle("Model 1a: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = veg.code.avg, y = shade.m * 100)) +
#   geom_line(data = df1a, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   xlab("Code for vegetation on both banks (0=grass; 1=forest)") +
#   ggtitle("Model 1a: Observed & Predicted Values")
# 
# #Linear regression model (Model #1b: Veg Type vs. Ln(Shade in Middle))
# mod1b <- lm(log(shade.m) ~ veg.code.avg, df.transect)
# anova(mod1b)  #Print ANOVA results
# summary(mod1b)  #Print coefficients
# df1b <- data.frame(x = vv, y = predict(mod1b, data.frame(veg.code.avg = vv)), 
#                    sres = rstandard(mod1b))  #df of standard residuals and predicted values
# ggplot(df1b, aes(x = x, y = sres)) + geom_point() +  #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Code for vegetation on both banks (0=grass; 1=forest)") + ylab("Standardized Residuals") +
#   ggtitle("Model 1b: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = veg.code.avg, y = log(shade.m))) +
#   geom_line(data = df1b, aes(x = x, y = y)) + ylab("ln(Shade in Middle of Stream)") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   xlab("Code for vegetation on both banks (0=grass; 1=forest)") +
#   ggtitle("Model 1b: Observed & Predicted Values")
# 
# #Multiple linear regression model (Model #2: all parameters for shade in middle of stream)
# mod2 <- lm(shade.m ~ wetted.width.m + thalweg.depth.m + lens.height.m + veg.code.avg + 
#              wetted.edge.l.m + wetted.edge.r.m + herb.height.max.avg.m + strm.azimuth + veg.height.max.l.m +
#              veg.height.max.r.m, df.transect)
# anova(mod2)  #Print ANOVA results
# summary(mod2)  #Print coefficients
# df2 <- data.frame(ypred = fitted(mod2), sres = rstandard(mod2)) 
# ggplot(df2, aes(x = ypred * 100, y = sres)) + geom_point() +  #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Predicted Shade in Middle of Transect (%)") + ylab("Standardized Residuals") +
#   ggtitle("Model 2: Standardized Residuals vs. Dependent Variable")
# 
# #Multiple linear regression model (Model #3: ind vars sig at 5% level for shade in middle of stream)
# mod3 <- lm(shade.m ~ wetted.width.m + thalweg.depth.m + lens.height.m + veg.code.avg + 
#              strm.azimuth + veg.dens.avg + herb.height.mode.avg.m, df.transect)
# anova(mod3)  #Print ANOVA results
# summary(mod3)  #Print coefficients
# df3 <- data.frame(ypred = fitted(mod3), sres = rstandard(mod3)) 
# ggplot(df3, aes(x = ypred * 100, y = sres)) + geom_point() +  #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Predicted Shade in Middle of Transect (%)") + ylab("Standardized Residuals") +
#   ggtitle("Model 3: Standardized Residuals vs. Dependent Variable")
# 
# #Multiple linear regression model (Model #4: ind vars sig at 5% level for shade in middle of stream)
# mod4 <- lm(shade.m ~ wetted.width.m + thalweg.depth.m + lens.height.m + veg.code.avg + 
#              veg.dens.avg, df.transect)
# anova(mod4)  #Print ANOVA results
# summary(mod4)  #Print coefficients
# df4 <- data.frame(ypred = fitted(mod4), sres = rstandard(mod4)) 
# ggplot(df4, aes(x = ypred * 100, y = sres)) + geom_point() +  #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Predicted Shade in Middle of Transect (%)") + ylab("Standardized Residuals") +
#   ggtitle("Model 4: Standardized Residuals vs. Dependent Variable")
# 
# #Linear regression model (Model #5: Wetted Width vs. Shade in Middle)
# mod5 <- lm(shade.m ~ wetted.width.m, df.transect)
# anova(mod5)  #Print ANOVA results
# summary(mod5)  #Print coefficients
# df5 <- data.frame(x = ww, y = predict(mod5, data.frame(wetted.width.m = ww)), 
#                   sres = rstandard(mod5))  #df of standard residual and predicted values from the linear model
# ggplot(df5, aes(x = x, y = sres)) + geom_point() +  #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Wetted Width (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 5: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = wetted.width.m, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df5, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("Wetted Width (m)") + ggtitle("Model 5: Observed & Predicted Values") +
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #6: Thalweg Depth vs. Shade in Middle)
# mod6 <- lm(shade.m ~ thalweg.depth.m, df.transect)
# anova(mod6)  #Print ANOVA results
# summary(mod6)  #Print coefficients
# df6 <- data.frame(x = tw, y = predict(mod6, data.frame(thalweg.depth.m = tw)), 
#                   sres = rstandard(mod6))  #df of predicted and standard residual values from the linear model
# ggplot(df6, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Thalweg Depth (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 6: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = thalweg.depth.m, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df6, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("Thalweg Depth (m)") + ggtitle("Model 6: Observed & Predicted Values") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #7: Lens Height Middle vs. Shade in Middle)
# mod7 <- lm(shade.m ~ lens.height.m, df.transect)
# anova(mod7)  #Print ANOVA results
# summary(mod7)  #Print coefficients
# df7 <- data.frame(x = lhm, y = predict(mod7, data.frame(lens.height.m = lhm)), 
#                   sres = rstandard(mod7))  #df of predicted and standard residual values from the linear model
# ggplot(df7, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Lens Height from Middle (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 7: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = lens.height.m, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df7, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("Lens Height from Middle (m)") + ggtitle("Model 7: Observed & Predicted Values") +
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #8: Lens Height Left vs. Shade in Left)
# mod8 <- lm(shade.l ~ lens.height.l, df.transect)
# anova(mod8)  #Print ANOVA results
# summary(mod8)  #Print coefficients
# df8 <- data.frame(x = lhl, y = predict(mod8, data.frame(lens.height.l = lhl)), 
#                        sres = rstandard(mod8))  #df of predicted and standard residual values from the linear model
# ggplot(df8, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Lens Height from Left (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 8: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = lens.height.l, y = shade.l * 100, colour = veg.code.avg)) +
#   geom_line(data = df8, aes(x = x, y = y * 100)) + ylab("% Shade on Left of Stream Transect") +
#   xlab("Lens Height from Left (m)") + ggtitle("Model 8: Observed & Predicted Values") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #9: Lens Height Right vs. Shade in Right)
# mod9 <- lm(shade.r ~ lens.height.r, df.transect)
# anova(mod9)  #Print ANOVA results
# summary(mod9)  #Print coefficients
# df9 <- data.frame(x = lhr, y = predict(mod9, data.frame(lens.height.r = lhr)),
#                   sres = rstandard(mod9))  #df of predicted and standard residual values from the linear model
# ggplot(df9, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Lens Height from Right (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 9: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = lens.height.r, y = shade.r * 100, colour = veg.code.avg)) +
#   geom_line(data = df9, aes(x = x, y = y * 100)) + ylab("% Shade on Right of Stream Transect") +
#   xlab("Lens Height from Right (m)") + ggtitle("Model 9: Observed & Predicted Values") +
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #10: Wetted Edge Left vs. Shade in Middle)
# mod10 <- lm(shade.m ~ wetted.edge.l.m, df.transect)
# anova(mod10)  #Print ANOVA results
# summary(mod10)  #Print coefficients
# df10 <- data.frame(x = wel, y = predict(mod10, data.frame(wetted.edge.l.m = wel)),
#                    sres = rstandard(mod10))  #create a dataframe with standard residual values from the linear model
# ggplot(df10, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Wetted Edge Left (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 10 (wetted edge left vs shade): Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = wetted.edge.l.m, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df10, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("Wetted Edge Left (m)") + ggtitle("Model 10: Observed & Predicted Values") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #11: Wetted Edge Right vs. Shade in Middle)
# mod11 <- lm(shade.m ~ wetted.edge.r.m, df.transect)
# anova(mod11)  #Print ANOVA results
# summary(mod11)  #Print coefficients
# df11 <- data.frame(x = wer, y = predict(mod11, data.frame(wetted.edge.r.m = wer)),
#                    sres = rstandard(mod11))  #create a dataframe with standard residual values from the linear model
# ggplot(df11, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Wetted Edge Left (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 11: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = wetted.edge.r.m, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df11, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("Wetted Edge Right (m)") + ggtitle("Model 11: Observed & Predicted Values") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #12: Wetted to depth vs. Shade in Middle)
# mod12 <- lm(shade.m ~ width.to.depth, df.transect)
# anova(mod12)  #Print ANOVA results
# summary(mod12)  #Print coefficients
# df12 <- data.frame(x = wwtw, y = predict(mod12, data.frame(width.to.depth = wwtw)), 
#                    sres = rstandard(mod12))  #df of predicted and standard residual values from the linear model
# ggplot(df12, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Width to Depth Ratio") + ylab("Standardized Residuals") +
#   ggtitle("Model 12: Standardized Residuals vs. Independent Variable")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = width.to.depth, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df12, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("Width to Depth Ratio") + ggtitle("Model 12: Observed & Predicted Values") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4))
# 
# #Linear regression model (Model #13: General Azimuth of Stream vs. Shade in Middle)
# mod13 <- lm(shade.m ~ gen.strm.azimuth, df.transect)
# anova(mod13)  #Print ANOVA results
# summary(mod13)  #Print coefficients
# df13 <- data.frame(x = azgen, y = predict(mod13, data.frame(gen.strm.azimuth = azgen)), 
#                    sres = rstandard(mod13))  #df of predicted and standard residual values from the linear model
# ggplot(df13, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("General Stream Azimuth (degrees)") + ylab("Standardized Residuals") +
#   ggtitle("Model 13: Standardized Residuals vs. Independent Variable") + 
#   labs(caption = "The departure angle of the stream from a south reference line (North-South = 0, Northwest-Southeast = -45)")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = gen.strm.azimuth, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df13, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("General Stream Azimuth (degrees)") + ggtitle("Model 13: Observed & Predicted Values") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4)) + 
#   labs(caption = "The departure angle of the stream from a south reference line (North-South = 0, Northwest-Southeast = -45)")
# 
# #Linear regression model (Model #14: Max Veg Height on Right vs. Shade in Middle)
# mod14 <- lm(shade.m ~ veg.height.max.r.m, df.transect)
# anova(mod14)  #Print ANOVA results
# summary(mod14)  #Print coefficients
# df14 <- data.frame(x = vhmaxr, y = predict(mod14, data.frame(veg.height.max.r.m = vhmaxr)), 
#                    sres = rstandard(mod14))  #df of predicted and standard residual values from the linear model
# ggplot(df14, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
#   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
#   xlab("Maximum Height of Vegetation on Right Bank (m)") + ylab("Standardized Residuals") +
#   ggtitle("Model 14: Standardized Residuals vs. Independent Variable") + 
#   labs(caption = "within a 10 m offset from the creek centerline")
# ggplot() +  #plot observed and predicted values
#   geom_point(data = df.transect, aes(x = veg.height.max.r.m, y = shade.m * 100, colour = veg.code.avg)) +
#   geom_line(data = df14, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
#   xlab("Maximum Height of Vegetation on Right Bank (m)") + ggtitle("Model 14: Observed & Predicted Values") + 
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   scale_colour_gradientn(colours = rainbow(4)) + 
#   labs(caption = "within a 10 m offset from the creek centerline")
# 
# 
