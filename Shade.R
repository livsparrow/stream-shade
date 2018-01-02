#Brown's Creek Riparian Shading Study December 22, 2017
#by Olivia Sparrow
#Civil Engineering Master's Student, University of Minnesota-Twin Cities
#Water Resources Engineer, Emmons & Olivier Resources

#Description of code: This code reads in the data I collected in Brown's Creek in 2017 and the outputs
#of processing hemispherical photographs using the software WinSCANOPY to estimate shade along the creek.


# Section 1: To Do List ---------------------------------------------------


#TO DO:
#label each transect on x axis
#add note on the assumed growing season
#import and plot shade over an entire season at one site
#fix variable labels in graphs to be more appropriate (i.e. not just their parameter numbers with periods)
#note that it's left to right looking downstream
#figure out how to import and use GIS data
#TO DO: Whisker box plots of variables?


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
library(wesanderson)  #another color option
library(classInt)
library(lidR)

# Section 3: Read and format data -----------------------------------------


#Read in and format the field data
df <- read.csv('TidyTransectData.csv', header = TRUE)  #Read full table of data collected
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
shp.boundary <- st_read("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/Boundary2.shp")  #Study Area Bounndary
shp.creek <- st_read("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/Creek.shp")  #Creek line
shp.shade <- st_read("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/BC_segs_shade_diss.shp")  #Thermal Study from Shade analysis
shp.shade.dem <- st_read("H:/2017 BCWD Riparian Shading Study/R/stream-shade/GIS/out_global_radiation2.shp")  #Topographic shade

#read in DEM shade analysis output and calculate topographic shade
shp.shade.dem$S4 <- 1 - shp.shade.dem$T4 / shp.shade.dem$T4[100]  #Monthly shade in May
shp.shade.dem$S5 <- 1 - shp.shade.dem$T5 / shp.shade.dem$T5[100]  #Monthly shade in June
shp.shade.dem$S6 <- 1 - shp.shade.dem$T6 / shp.shade.dem$T6[100]  #Monthly shade in July
shp.shade.dem$S7 <- 1 - shp.shade.dem$T7 / shp.shade.dem$T7[100]  #Monthly shade in August
shp.shade.dem$S8 <- 1 - shp.shade.dem$T8 / shp.shade.dem$T8[100]  #Monthly shade in September
shp.shade.dem$growing.season.shade <- (shp.shade.dem$S4 + shp.shade.dem$S5 +  #Average shade over growing season
  shp.shade.dem$S6 + shp.shade.dem$S7 + shp.shade.dem$S8) / 5
ss.dem <- as.data.frame(shp.shade.dem$growing.season.shade[1:99])  #convert topographic shade to data frame
df$shade.lidar.dem <- ss.dem  #save topographic shade to df

#TO DO: Use LiDAR processing when finished in Arcmap or separate code
#save to df at transect points
#may need to access heigh along entire stream for modeling

#add northings and eastings to df since coordinates are only as lat and long format
utm15nCRS <- st_crs(shp.shade)  #save geospatial metadata
cord.dec <- SpatialPoints(cbind(df$long.dec.deg, df$lat.dec.deg), 
                          proj4string = CRS("+proj=longlat"))  #save existing coodinates in lat-long system
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:26915"))
df.cord.utm <- as.data.frame(cord.UTM)
df$easting <- df.cord.utm$coords.x1
df$northing <- df.cord.utm$coords.x2
write.csv(df, "data_table_with_coords") #write df to csv for use in ArcGIS





# Section 4: Spatial Data Analysis --------------------------------------------

#plot the project boundary and creek location (Note that coordinates are in Northing and Easting)
plot(st_geometry(shp.boundary), border = "black", lwd = 2, main = "Study Area Boundary")  #plots boundary of study area
plot(st_geometry(shp.creek), add = TRUE, col = "blue")

#plot the shade estimated by GIS analysis of 2011 LiDAR (Note that coordinates are in Northing and Easting)
pal2 <- rev(wes_palette(n = 5, name = "Zissou"))
br <- c(0, 20, 40, 60, 80, 100)
offs <- 0.0000001 
br[1] <- br[1] - offs 
br[length(br)] <- br[length(br)] + offs 
shade_cuts <- cut(shp.shade$Shade, br)
plot(shp.shade["Shade"], lwd = 2, axes = TRUE, col = pal2[as.numeric(shade_cuts)],
     xlab = "Easting", ylab = "Northing", 
     main = "Shade from GIS Analysis of 2011 LiDAR ")  #plot shade estimated using ArcGIS Solar Radiation tool and LiDAR
legend("topright", legend = paste("<", round(br[-1])), col = pal2, lty = 1, lwd = 2)

#add transect locations to the above plot (Note that coordinates are in lat and lon)
pts.transects <- st_as_sf(df, coords = c("easting", "northing"), crs = utm15nCRS)  #convert field data to sf for map
plot(shp.shade["Shade"], lwd = 2, axes = TRUE, col = pal2[as.numeric(shade_cuts)],
     xlab = "Easting", ylab = "Northing", 
     main = "Shade from GIS Analysis of 2011 LiDAR & Riparian Shading Study Transects")  #plot shade estimated using ArcGIS Solar Radiation tool and LiDAR
legend("topright", legend = c("Transects", paste("<", round(br[-1]))), pch = c(1, NA, NA, NA, NA, NA, NA), 
       lty = c(NA, 1, 1, 1, 1, 1, 1), col = c("black", pal2), lwd = 2)
plot(st_geometry(pts.transects), add = TRUE)  # add transect locations to plot

#determine and save the shade estimated using lidar into the df
dist <- st_distance(pts.transects, shp.shade)  #calculate distance between each creek segment and transect point
i <- as.matrix(apply(dist, 1, which.min))  #find the row index in shp.shade.buf with the shortest distance to each point
ss <- shp.shade$Shade[i]
pts.transects$shade.lidar.leafoff.dsm <- ss
df[, "shade.lidar.leafoff.dsm"] <- ss
shade_cuts <- cut(pts.transects$shade.lidar.leafoff.dsm, br)
plot(pts.transects["shade.lidar.leafoff.dsm"], axes = TRUE, col = pal2[as.numeric(shade_cuts)],
     xlab = "Easting", ylab = "Northing", 
     main = "Shade from GIS Analysis of 2011 LiDAR at Riparian Shading Study Transects")
legend("topright", legend = paste("<", round(br[-1])), col = pal2, pch = 1, lwd = 2)







#shp.shade.lidar <- readOGR("H:/2017 BCWD Riparian Shading Study/R/Shade/GIS/BC_Segs_Shade.shp", "BC_Segs_Shade")
#spplot(shp.creek)
#and use to populate df

# Section 5: Subset Dataframes -------------------------------------------


#subset the data into tables for use in analysis
df.transect <- subset(df, purpose == 'Transect')  #Create a smaller table of just the transect observations
df.transect.num <- df.transect[, sapply(df.transect, is.numeric)]  #Create df of only numeric values for correlation analysis
df.transect.grass <- subset(df.transect, veg.code.avg == 0)  #Numerical observations at grassy & mixed transects
df.stage <- subset(df, purpose == 'Stage')  #Create a smaller table of just the staging observations
df.test <- subset(df,purpose == 'Test')  #Create a smaller table of just the two sites to be used to test the regression
#str(df.transect)  #Print structure of data collection table
#str(df.transect.grass)  #Print structure of data collection table

#Create long dataframes for plotting results for each transect and stage-shade curves
df.transect.m <- melt(df.transect, id.vars = c("reach.id", "transect.no"), 
                   measure.vars = c("shade.l", "shade.m", "shade.r", "shade.avg"))
#Break the stage and height dataframes into separate objects so they can be melted and then merged later
df.stage.height <- data.frame("reach.id" = df.stage$reach.id, "transect.no" = df.stage$transect.no,
                              "Left" = df.stage$lens.height.l, "Middle" = df.stage$lens.height.m,  
                              "Right" = df.stage$lens.height.r)
df.stage.shade <- data.frame("reach.id" = df.stage$reach.id, "transect.no" = df.stage$transect.no,
                              "Left" = df.stage$shade.l, "Middle" = df.stage$shade.m,  
                              "Right" = df.stage$shade.r)
#melt and rename stage (height) and shade data frames by position
df.stage.height.m <- melt(df.stage.height, id.vars = c("reach.id", "transect.no"), 
                          measure.vars = c("Left", "Middle", "Right"))
colnames(df.stage.height.m) <- c("reach.id", "transect.no", "position", "lens.height")
df.stage.shade.m <- melt(df.stage.shade, id.vars = c("reach.id", "transect.no"), 
                          measure.vars = c("Left", "Middle", "Right"))
colnames(df.stage.shade.m) <- c("reach.id", "transect.no", "position", "shade")
#merge staging dataframes
df.stage.m <- data.frame(df.stage.height.m, "shade" = df.stage.shade.m$shade)
setorder(df.stage.m, lens.height)
head(df.stage.m)
#subset staging df by transect and position
df.stage.m.36L <- subset(df.stage.m, transect.no == 6 & position == "Left")
df.stage.m.36M <- subset(df.stage.m, transect.no == 6 & position == "Middle")
df.stage.m.36R <- subset(df.stage.m, transect.no == 6 & position == "Right")
df.stage.m.39L <- subset(df.stage.m, transect.no == 9 & position == "Left")
df.stage.m.39M <- subset(df.stage.m, transect.no == 9 & position == "Middle")
df.stage.m.39R <- subset(df.stage.m, transect.no == 9 & position == "Right")
#write each to curve CSV
write.csv(df.stage.m, "staging tables")


# Section 6: Regression of Stage-Shade Curves and Correction --------------


#Plot shade-stage curves for 4 reference transects
ggplot(df.stage.m, aes(lens.height, shade * 100, colour = position)) + geom_point() + geom_line() +
  labs(title = "1a: Stage-Shade Curves from WinSCANOPY Simulation", x = "Height of lens above water (m)", 
       y = "Average % Shade Over Growing Season") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  facet_wrap(~ transect.no,  labeller = label_both)
      #TODO: label reaches in each grid
      #TODO: check why bottom right reach has such as large difference in left position shade

#METHOD 1: Using standard x value, interpolate corresponding y value
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
  df.ys[i, 4] <- min(df.stage.m$lens.height[df.stage.m$reach.id == r & df.stage.m$transect.no == t &
                                     df.stage.m$position == p])
  df.ys[i, 5] <- min(df.stage.m$lens.height[df.stage.m$reach.id == r & df.stage.m$transect.no == t &
                                     df.stage.m$position == p & df.stage.m$lens.height > df.ys[i, 4]])
  df.ys[i, 6] <- df.stage.m$shade[df.stage.m$reach.id == r & df.stage.m$transect.no == t &
                                     df.stage.m$position == p & df.stage.m$lens.height == df.ys[i, 4]]
  df.ys[i, 7] <- df.stage.m$shade[df.stage.m$reach.id == r & df.stage.m$transect.no == t &
                           df.stage.m$position == p & df.stage.m$lens.height == df.ys[i, 5]]

}
#interpolate reference value for shade
df.ys$ys <- (df.ys$xs - df.ys$x1) * (df.ys$y2 - df.ys$y1) / (df.ys$x2 - df.ys$x1) + df.ys$y1  
names(df.ys)[1] <- "reach.id"
names(df.ys)[2] <- "transect.no"
names(df.ys)[3] <- "position"
#save reference values to staging df
df.stage.mer <- merge(df.stage.m, df.ys, by = c("reach.id", "transect.no","position"), all.x = TRUE)
#Standardize the height and shade using reference values
df.stage.mer$xstar1 <- df.stage.mer$lens.height / df.stage.mer$xs  #standardize lens height above stream
df.stage.mer$ystar1 <- df.stage.mer$shade / df.stage.mer$ys  #standardize shade
#Plot normalized shade-stage curves for 2 grassy reference transects
ggplot(df.stage.mer, aes(xstar1, ystar1, colour = position)) + geom_point() + geom_line() +
  labs(title = "1b: Standardized Stage-Shade Curves (WinSCANOPY) - Method 1", x = "Normalized Lens Height Above Water",
       y = "Normalized Average Shade Over Growing Season")+
  facet_wrap(~ transect.no,  labeller = label_both)

#Method 2: Standardize the height and shade values using standard normal deviate
#Loop through rows of stage dataframe to calculate standard normal deviate
c <- ncol(df.stage.mer)  #number of columns in df.stage.mer
for(i in 1:nrow(df.stage.mer)) #each i is a row in the input table.each measurement has its own rowthere are multiple rows and measurements in each transect
{
  #if(df.stage.mer[i, 1] == 3) {
  t <- df.stage.mer[i, 2]  #current transect
  p <- df.stage.mer[i, 3]  #current position (left, middle, right)
  #lens height:
  xbar <- mean(df.stage.mer$lens.height[df.stage.mer$transect.no == t & df.stage.mer$position == p])   #average lens height
  xs <- sd(df.stage.mer$lens.height[df.stage.mer$transect.no == t & df.stage.mer$position == p])    #standard deviation of height
  df.stage.mer[i, c + 1] <- (df.stage.mer[i, 4] - xbar) / xs   #standard normal deviate of measurement i
  #shade:
  ybar <- mean(df.stage.mer$shade[df.stage.mer$transect.no == t & df.stage.mer$position == p])   #average shade
  ys <- sd(df.stage.mer$shade[df.stage.mer$transect.no == t & df.stage.mer$position == p])    #standard deviation of shade
  df.stage.mer[i, c + 2] <- (df.stage.mer[i, 5] - ybar) / ys   #standard normal deviate of estimate i
  #}
}
names(df.stage.mer)[c + 1] <- "xstar2"
names(df.stage.mer)[c + 2] <- "ystar2"
#Plot normalized shade-stage curves for 4 reference transects
ggplot(df.stage.mer, aes(xstar2, ystar2, colour = position)) + geom_point() + geom_line() +
  labs(title = "1c: Standardized Stage-Shade Curves (WinSCANOPY) - Method 2", x = "Std norm dev of lens height above water", 
       y = "Std norm dev of Average % Shade Over Growing Season")+
  facet_wrap(~ transect.no,  labeller = label_both)

#TODO: once normalization method is finished - apply to df in order to correct estimated shade based on height above stream


# Section 7: Correlation Matrix -------------------------------------------


#correlation between independent & dependent variables
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

#TODO: Update these to plot updated shade based on height above stream
#Plot shade for each reach and transect and position
ggplot(df.transect.m, aes(transect.no, 100 * value, colour = variable)) + geom_point() + geom_line() +
  labs(title = "1. Estimated Shade from WinSCANOPY Simulation", x = "Transect #", 
       y = "Average % Shade Over Growing Season") +
  scale_x_continuous(breaks = seq(0, 11, 1)) + scale_y_continuous(breaks = seq(0, 100, 20)) +
  facet_wrap(~ reach.id,  labeller = label_both)

#Plot average shade at each transect and reach
ggplot(df.transect, aes(transect.no, shade.avg * 100, group = factor(reach.id), colour = factor(reach.id))) +
  geom_point() + geom_line() + scale_x_continuous(breaks = seq(0, 11, 1)) + scale_y_continuous(breaks = seq(0, 100, 20)) +
  labs(title = "2a. Average Shade at Each Transect from WinSCANOPY Simulation",
       x = "Transect #", y = "Average % Shade Over Growing Season") +
  facet_wrap(~ reach.id,  labeller = label_both)

#Plot average shade at each transect and reach - comparing hemiphotos to LiDAR
df.transect$shade.lidar.leafoff.dsm <- df.transect$shade.lidar.leafoff.dsm / 100
df.transect.shade.m <- melt(df.transect, id.vars = c("reach.id", "transect.no"), 
                      measure.vars = c("shade.avg", "shade.lidar.leafoff.dsm"))
colnames(df.transect.shade.m)[3] <- "method"  #update to reflect that variable column is the method of estimating shade
levels(df.transect.shade.m$method) <- c(levels(df.transect.shade.m$method), "WinSCANOPY", "LiDAR")
df.transect.shade.m$method[df.transect.shade.m$method == "shade.avg"] <- "WinSCANOPY"
df.transect.shade.m$method[df.transect.shade.m$method == "shade.lidar.leafoff.dsm"] <- "LiDAR"
ggplot(df.transect.shade.m, aes(transect.no, 100 * value, group = factor(reach.id), colour = method)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 11, 1)) + scale_y_continuous(breaks = seq(0, 100, 20)) +
  labs(title = "2b. Average Shade at Each Transect from WinSCANOPY Simulation and LiDAR Analysis",
       x = "Transect #", y = "Average % Shade Over Growing Season") +
  facet_wrap(~ reach.id,  labeller = label_both)


#Plot observed values vs. shade in middle of stream for grassy transects
df.transect.num.grass.m <- melt(df.transect.num.grass, id = c('shade.m'))
ggplot(data = df.transect.num.grass.m, aes(x = value, y = shade.m * 100, colour = variable)) +
  geom_point() +  xlab("Independent Variables") + ylab("% Shade in middle of stream transect") +
  ggtitle("6a. Observed Values (Grassy Transects)") + 
  facet_wrap( ~ variable, scales="free_x") + theme(legend.position = "none")

#Plot observed values vs. shade in middle of stream for grassy transects
df.transect.num.m <- melt(df.transect.num, id.vars = c('shade.m'))
ggplot(data = df.transect.num.m, aes(x = value, y = shade.m * 100, colour = variable)) +
  geom_point() +  xlab("Independent Variables") + ylab("% Shade in middle of stream transect") +
  ggtitle("6b. Observed Values (All Transects)") + scale_y_continuous(breaks = seq(0, 100, 20)) +
  facet_wrap( ~ variable, scales="free_x") + theme(legend.position = "none")

#Plot vegetation height for all transects wrapped by veg code avg
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

#Plot lens height based on vegetation type
df.transect.lens.m <- melt(df.transect, id.vars = c('reach.id', 'transect.no', 'veg.code.avg'), 
                           measure.vars = c('lens.height.m', 'lens.height.l', 'lens.height.r'))
df.transect.lens.m$veg.code.avg <- factor(df.transect.lens.m$veg.code.avg)
ggplot(df.transect.lens.m, aes(veg.code.avg, value)) + geom_boxplot() + xlab("Riparian Vegetation (Grass = 0...Forest = 1)") +
  ylab("Lens Height Above Water (m)") + ggtitle("8. Height of Camera Lens Above Stream (All Transects)")

# Section 9: Regression Models --------------------------------------------

#Create vectors for the independent variables
vv <- df.transect$veg.code.avg
ww <- df.transect$wetted.width.m 
tw <- df.transect$thalweg.depth.m
wwtw <- df.transect$width.to.depth
lhm <- df.transect$lens.height.m
lhl <- df.transect$lens.height.l  
lhr <- df.transect$lens.height.r
wel <- df.transect$wetted.edge.l.m
wer <- df.transect$wetted.edge.r.m
hmaxavg <- df.transect$herb.height.max.avg.m 
asp <- df.transect$stream.aspect.deg
ori <- df.transect$stream.orient.from.north
vdens <- df.transect$veg.dens.avg
hmodavg <- df.transect$herb.height.mode.avg.m
#ypred = predict(mod2, data.frame(wetted.width.m = ww, thalweg.depth.m = tw, lens.height.m = lhm, 
#veg.code.avg = vv, wetted.edge.l.m = wel, wetted.edge.r.m = wer, 
#herb.height.max.avg.m = hmaxavg, stream.aspect.deg = asp, 
#veg.dens.avg = vdens, herb.height.mode.avg.m = hmodavg)))

#Linear regression model (Model #1a: Veg Type vs. Shade in Middle)
mod1a <- lm(shade.m ~ veg.code.avg, df.transect)
anova(mod1a)  #Print ANOVA results
summary(mod1a)  #Print coefficients
df1a <- data.frame(x = vv, y = predict(mod1a, data.frame(veg.code.avg = vv)), 
                   sres = rstandard(mod1a))  #df of standard residuals and predicted values
ggplot(df1a, aes(x = x, y = sres)) + geom_point() +  #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Code for vegetation on both banks (0=grass; 1=forest)") + ylab("Standardized Residuals") +
  ggtitle("Model 1a: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = veg.code.avg, y = shade.m * 100)) +
  geom_line(data = df1a, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  xlab("Code for vegetation on both banks (0=grass; 1=forest)") +
  ggtitle("Model 1a: Observed & Predicted Values")

#Linear regression model (Model #1b: Veg Type vs. Ln(Shade in Middle))
mod1b <- lm(log(shade.m) ~ veg.code.avg, df.transect)
anova(mod1b)  #Print ANOVA results
summary(mod1b)  #Print coefficients
df1b <- data.frame(x = vv, y = predict(mod1b, data.frame(veg.code.avg = vv)), 
                   sres = rstandard(mod1b))  #df of standard residuals and predicted values
ggplot(df1b, aes(x = x, y = sres)) + geom_point() +  #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Code for vegetation on both banks (0=grass; 1=forest)") + ylab("Standardized Residuals") +
  ggtitle("Model 1b: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = veg.code.avg, y = log(shade.m))) +
  geom_line(data = df1b, aes(x = x, y = y)) + ylab("ln(Shade in Middle of Stream)") + 
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  xlab("Code for vegetation on both banks (0=grass; 1=forest)") +
  ggtitle("Model 1b: Observed & Predicted Values")

#Multiple linear regression model (Model #2: all parameters for shade in middle of stream)
mod2 <- lm(shade.m ~ wetted.width.m + thalweg.depth.m + lens.height.m + veg.code.avg + 
             wetted.edge.l.m + wetted.edge.r.m + herb.height.max.avg.m + stream.aspect.deg, df.transect)
anova(mod2)  #Print ANOVA results
summary(mod2)  #Print coefficients
df2 <- data.frame(ypred = fitted(mod2), sres = rstandard(mod2)) 
ggplot(df2, aes(x = ypred * 100, y = sres)) + geom_point() +  #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Predicted Shade in Middle of Transect (%)") + ylab("Standardized Residuals") +
  ggtitle("Model 2: Standardized Residuals vs. Dependent Variable")

#Multiple linear regression model (Model #3: ind vars sig at 5% level for shade in middle of stream)
mod3 <- lm(shade.m ~ wetted.width.m + thalweg.depth.m + lens.height.m + veg.code.avg + 
            stream.aspect.deg + veg.dens.avg + herb.height.mode.avg.m, df.transect)
anova(mod3)  #Print ANOVA results
summary(mod3)  #Print coefficients
df3 <- data.frame(ypred = fitted(mod3), sres = rstandard(mod3)) 
ggplot(df3, aes(x = ypred * 100, y = sres)) + geom_point() +  #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Predicted Shade in Middle of Transect (%)") + ylab("Standardized Residuals") +
  ggtitle("Model 3: Standardized Residuals vs. Dependent Variable")

#Multiple linear regression model (Model #4: ind vars sig at 5% level for shade in middle of stream)
mod4 <- lm(shade.m ~ wetted.width.m + thalweg.depth.m + lens.height.m + veg.code.avg + 
             veg.dens.avg, df.transect)
anova(mod4)  #Print ANOVA results
summary(mod4)  #Print coefficients
df4 <- data.frame(ypred = fitted(mod4), sres = rstandard(mod4)) 
ggplot(df4, aes(x = ypred * 100, y = sres)) + geom_point() +  #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Predicted Shade in Middle of Transect (%)") + ylab("Standardized Residuals") +
  ggtitle("Model 4: Standardized Residuals vs. Dependent Variable")

#Linear regression model (Model #5: Wetted Width vs. Shade in Middle)
mod5 <- lm(shade.m ~ wetted.width.m, df.transect)
anova(mod5)  #Print ANOVA results
summary(mod5)  #Print coefficients
df5 <- data.frame(x = ww, y = predict(mod5, data.frame(wetted.width.m = ww)), 
                  sres = rstandard(mod5))  #df of standard residual and predicted values from the linear model
ggplot(df5, aes(x = x, y = sres)) + geom_point() +  #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Wetted Width (m)") + ylab("Standardized Residuals") +
  ggtitle("Model 5: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = wetted.width.m, y = shade.m * 100, colour = veg.code.avg)) +
  geom_line(data = df5, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  xlab("Wetted Width (m)") + ggtitle("Model 5: Observed & Predicted Values") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #6: Thalweg Depth vs. Shade in Middle)
mod6 <- lm(shade.m ~ thalweg.depth.m, df.transect)
anova(mod6)  #Print ANOVA results
summary(mod6)  #Print coefficients
df6 <- data.frame(x = tw, y = predict(mod6, data.frame(thalweg.depth.m = tw)), 
                  sres = rstandard(mod6))  #df of predicted and standard residual values from the linear model
ggplot(df6, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Thalweg Depth (m)") + ylab("Standardized Residuals") +
  ggtitle("Model 6: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = thalweg.depth.m, y = shade.m * 100, colour = veg.code.avg)) +
  geom_line(data = df6, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  xlab("Thalweg Depth (m)") + ggtitle("Model 6: Observed & Predicted Values") + 
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #7: Lens Height Middle vs. Shade in Middle)
mod7 <- lm(shade.m ~ lens.height.m, df.transect)
anova(mod7)  #Print ANOVA results
summary(mod7)  #Print coefficients
df7 <- data.frame(x = lhm, y = predict(mod7, data.frame(lens.height.m = lhm)), 
                  sres = rstandard(mod7))  #df of predicted and standard residual values from the linear model
ggplot(df7, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Lens Height from Middle (m)") + ylab("Standardized Residuals") +
  ggtitle("Model 7: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = lens.height.m, y = shade.m * 100, colour = veg.code.avg)) +
  geom_line(data = df7, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  xlab("Lens Height from Middle (m)") + ggtitle("Model 7: Observed & Predicted Values") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #8: Lens Height Left vs. Shade in Left)
mod8 <- lm(shade.l ~ lens.height.l, df.transect)
anova(mod8)  #Print ANOVA results
summary(mod8)  #Print coefficients
df8 <- data.frame(x = lhl, y = predict(mod8, data.frame(lens.height.l = lhl)), 
                       sres = rstandard(mod8))  #df of predicted and standard residual values from the linear model
ggplot(df8, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Lens Height from Left (m)") + ylab("Standardized Residuals") +
  ggtitle("Model 8: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = lens.height.l, y = shade.l * 100, colour = veg.code.avg)) +
  geom_line(data = df8, aes(x = x, y = y * 100)) + ylab("% Shade on Left of Stream Transect") +
  xlab("Lens Height from Left (m)") + ggtitle("Model 8: Observed & Predicted Values") + 
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #9: Lens Height Right vs. Shade in Right)
mod9 <- lm(shade.r ~ lens.height.r, df.transect)
anova(mod9)  #Print ANOVA results
summary(mod9)  #Print coefficients
df9 <- data.frame(x = lhr, y = predict(mod9, data.frame(lens.height.r = lhr)),
                  sres = rstandard(mod9))  #df of predicted and standard residual values from the linear model
ggplot(df9, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Lens Height from Right (m)") + ylab("Standardized Residuals") +
  ggtitle("Model 9: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = lens.height.r, y = shade.r * 100, colour = veg.code.avg)) +
  geom_line(data = df9, aes(x = x, y = y * 100)) + ylab("% Shade on Right of Stream Transect") +
  xlab("Lens Height from Right (m)") + ggtitle("Model 9: Observed & Predicted Values") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #10: Wetted Edge Left vs. Shade in Middle)
mod10 <- lm(shade.m ~ wetted.edge.l.m, df.transect)
anova(mod10)  #Print ANOVA results
summary(mod10)  #Print coefficients
df10 <- data.frame(x = wel, y = predict(mod10, data.frame(wetted.edge.l.m = wel)),
                   sres = rstandard(mod10))  #create a dataframe with standard residual values from the linear model
ggplot(df10, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Wetted Edge Left (m)") + ylab("Standardized Residuals") +
  ggtitle("Model 10 (wetted edge left vs shade): Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = wetted.edge.l.m, y = shade.m * 100, colour = veg.code.avg)) +
  geom_line(data = df10, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  xlab("Wetted Edge Left (m)") + ggtitle("Model 10: Observed & Predicted Values") + 
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #11: Wetted Edge Right vs. Shade in Middle)
mod11 <- lm(shade.m ~ wetted.edge.r.m, df.transect)
anova(mod11)  #Print ANOVA results
summary(mod11)  #Print coefficients
df11 <- data.frame(x = wer, y = predict(mod11, data.frame(wetted.edge.r.m = wer)),
                   sres = rstandard(mod11))  #create a dataframe with standard residual values from the linear model
ggplot(df11, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Wetted Edge Left (m)") + ylab("Standardized Residuals") +
  ggtitle("Model 11: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = wetted.edge.r.m, y = shade.m * 100, colour = veg.code.avg)) +
  geom_line(data = df11, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  xlab("Wetted Edge Right (m)") + ggtitle("Model 11: Observed & Predicted Values") + 
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #12: Wetted to depth vs. Shade in Middle)
mod12 <- lm(shade.m ~ width.to.depth, df.transect)
anova(mod12)  #Print ANOVA results
summary(mod12)  #Print coefficients
df12 <- data.frame(x = wwtw, y = predict(mod12, data.frame(width.to.depth = wwtw)), 
                   sres = rstandard(mod12))  #df of predicted and standard residual values from the linear model
ggplot(df12, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Width to Depth Ratio") + ylab("Standardized Residuals") +
  ggtitle("Model 12: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = width.to.depth, y = shade.m * 100, colour = veg.code.avg)) +
  geom_line(data = df12, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  xlab("Width to Depth Ratio") + ggtitle("Model 12: Observed & Predicted Values") + 
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))

#Linear regression model (Model #13: Orientation of Stream vs. Shade in Middle)
mod13 <- lm(shade.m ~ stream.orient.from.north, df.transect)
anova(mod13)  #Print ANOVA results
summary(mod13)  #Print coefficients
df13 <- data.frame(x = ori, y = predict(mod13, data.frame(stream.orient.from.north = ori)), 
                   sres = rstandard(mod13))  #df of predicted and standard residual values from the linear model
ggplot(df13, aes(x = x, y = sres)) + geom_point() + #Standardized residuals plot
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Orientation of Stream (degrees)") + ylab("Standardized Residuals") +
  ggtitle("Model 13: Standardized Residuals vs. Independent Variable")
ggplot() +  #plot observed and predicted values
  geom_point(data = df.transect, aes(x = stream.orient.from.north, y = shade.m * 100, colour = veg.code.avg)) +
  geom_line(data = df13, aes(x = x, y = y * 100)) + ylab("% Shade in Middle of Stream") +
  xlab("Orientation of Stream (degrees)") + ggtitle("Model 13: Observed & Predicted Values") + 
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_gradientn(colours = rainbow(4))



# Section 10: Test Regression Model ---------------------------------------


#Read short table of only the transects that were repeated (these are the first site visits - the second were in the stage-storage curves)
#shaderep=read.csv('WinSCANOPYtransect_shade_repeat.csv',header=TRUE)

#TO DO: Fix the plots below
#Plot repeated visits on shade-stage curves for 4 reference transects
#ggplot()+geom_point(data=stage,aes(Height.m,Shade*100,colour=Position))+
#  geom_point(data=shaderep,aes(Height.m,Shade*100,colour=Position))+
#  labs(title="Stage-Shade Curves from WinSCANOPY Simulation",x="Height of lens above water (m)",y = "% Shade")+
#  facet_wrap(~Transect,  labeller = label_both)
