###############################################
#       Pilbara Toad Density Mapping          #
###############################################

#### Clear environment and set working directory
rm(list = ls())

# load libraries
library(sf)
library(raster)
library(sp)
library(dplyr)
library(fields)
library(gstat)


# Add current date to environment for output of file names with date
currentDate <- Sys.Date()

#prevent display in scientific notation
options(scipen = 999)

# Read in all perennial water point datasets
# These were previously filtered by Brenton in QGIS, originating from layers that Judy provided

dsn <- "."
lay <- "waterbody_pointsFRESH"
bodypoints <- read_sf(dsn = dsn, layer = lay)
bodypoints <- st_cast(bodypoints, "MULTIPOINT")

lay <- "clippedwater_points"
waterpoints <- read_sf(dsn = dsn, layer = lay)
waterpoints <- st_cast(waterpoints, "MULTIPOINT")

lay <- "clip_and_filt_drainage_points"
waterholes <- read_sf(dsn = dsn, layer = lay)
waterholes <- st_cast(waterholes, "MULTIPOINT")

# Combine into single object
bodypoints <- coordinates(as_Spatial(bodypoints))
waterpoints <- coordinates(as_Spatial(waterpoints))
waterholes <- coordinates(as_Spatial(waterholes))

allpoints <- as.data.frame(rbind(bodypoints,waterpoints,waterholes))
rm(bodypoints,waterpoints,waterholes)
colnames(allpoints) <- c("longitude","latitude")
head(allpoints)
coordinates(allpoints) <- ~longitude+latitude
plot(allpoints)

# Read in the number of rainy days (> 1 mm) for Australia (BOM dataset)
raindays <- raster("rain1mman/rdann-1.asc")

# extract the values of rainy days for each water point in the Pilbara
pointraindays <- extract(raindays, allpoints, sp=TRUE)
df <- as.data.frame(pointraindays)
colnames(df) <- c("Raindays","longitude","latitude")
hist(df$Raindays)
head(df)


## Adjust these Raindays values as per Tingley et al. 2013
# i.e. ndays=x+3(x)(1-(3(d-d^2)+d^3)), where d = (x-1)/364
df$Raindays <- df$Raindays+3*(df$Raindays)*(1-(3*(((df$Raindays-1)/364)-((df$Raindays-1)/364)^2)+((df$Raindays-1)/364)^3))
df$Raindays <- round(df$Raindays)

# Ignore the Tingley adjustment for big rainfall years as we're looking at the average


## Get the kernels for each combination of location and number of rainy days
# Load the kernel data
load("data/Kernel_fits.RData")
fits <- as.data.frame(fits)
head(fits)

# match the # of rainy days in the Pilbara point df with the row number in the kernel data
fits$Raindays <- as.numeric(rownames(fits))
df <- merge(df, fits, by = "Raindays")
head(df)


# For each point, j:
# 1. get pairwise distances to all water point
# 2. then calculate the density at j accruing from each of the water points

# # Function for the 2D kernel (x should be distance in metres)
 dcncross2D <- function(x, u, v) {
   (x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v)))/(2*pi*x)
 }


# Average perimeter length of 663 perennial waterbodies in the Pilbara is 1886.649 m
# Smart et al. 2020 calculated 350 adult toads per km = 0.350 per m
# Meaning the average water feature has 0.35*1886.649 = 660.3271 toads
# Probability density values then get multiplied by 660.3271 to get density of toads per m2
# Then this number gets multiplied by 1,000,000 to get number of toads per square km
# This number also aligns nicely with Florance et al's number of toads around AWPS: 672

# get coordinates of all waterpoints
all_coords <- df %>%
  dplyr::select(longitude, latitude) %>%
  as.matrix()


# Create SpatialPointsDataFrame from the toad density data
# toaddens <- grid.points #df[c(2,3,7)]
# coordinates(toaddens) <- c("longitude", "latitude")
# # Set the coordinate reference system (CRS) - replace EPSG:4326 with your CRS if different
# proj4string(toaddens) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Get another SpatialPointsDataFrame in which to dump the density values
# First sort out a Pilbara raster
austemp <- raster("annualmeantemp.tif")

lay <- "Pilbara_boundary_only"
pilbara <- read_sf(dsn = dsn, layer = lay)
# clip raster to Pilbara outline
pilbtemp <- crop(austemp, extent(pilbara))
pilbtemp <- mask(pilbtemp, pilbara)
# check it worked
plot(pilbtemp)
# turn this into a useable surface for interpolation
#pilbtempsp <- as(pilbtemp, 'SpatialPointsDataFrame')
pilbtempdf <- data.frame(as.data.frame(pilbtemp), coordinates(pilbtemp))


# Create a function to calculate the density at a given grid point 'j' from all water points
# first a vector to hold results
pointDensity <- rep(NA, nrow(pilbtempdf))

# run loop through each row and calculate density at that point
for (j in 1:nrow(pilbtempdf)) {
  if (is.na(pilbtempdf[j, 1])) next
  # Get the coordinates of the focal point 'j'
  focal_coords <- unlist(pilbtempdf[j, 2:3]) #c(grid.points$longitude[j], grid.points$latitude[j])
  
  
  # Calculate the distances in km from the focal grid point to all water points
  distances <- spDistsN1(all_coords, focal_coords, longlat = TRUE)
  
  # Convert distances to metres
  distances <- distances*1000
  
  # get the density at j from each water point
  dens <- dcncross2D(distances, df$u, df$v)
  
  # add sum of densities to correct column and row of df
  pointDensity[j] <- sum(dens, na.rm=TRUE)

}

# Convert these density values to density of toads per m2 (based on mean perimeter length of all waterpoints)
pointDensity <- pointDensity*660.3271 # Multiply density by number of toads at an average waterpoint

# Convert these density values to toads per km^2
pointDensity <- pointDensity*1000000

# Check that we've got some varying densities in our data frame

hist(pointDensity, breaks=100) #, ylim=c(0,6000)

# sanity check - have we got the same number of toads (approximately) as we started with?
sum(pointDensity, na.rm =TRUE) # given 1km^2 grid cells = total toads
660.3271 * nrow(all_coords) # number we started with
# pretty good: we lost a few to the ocean.



# replace raster values with calculated densities
pilbtemp@data@values <- pointDensity

plot(pilbtemp)

# show places where toad density > 1 per square kilometre
temp <- pilbtemp
temp@data@values <- as.numeric(temp@data@values>1)
plot(temp)



## Output per square km raster as a geotiff
writeRaster(pilbtemp, paste0("figures/toad_density_raster_km2_counts_", currentDate, ".tif"), format="GTiff", overwrite=TRUE)
writeRaster(temp, paste0("figures/toad_presence_raster_km2_counts_", currentDate, ".tif"), format="GTiff", overwrite=TRUE)
# writeRaster(toad_dens_raster, "GIS/BrentonStuff/toad_density_raster_20231031.tif", format="GTiff", overwrite=TRUE)



