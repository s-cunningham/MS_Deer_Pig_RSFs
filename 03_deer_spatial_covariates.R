library(tidyverse)
library(terra)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
deer <- read_csv("output/deer_used_avail_locations.csv") 

# Rasters
rast_list <- c("data/landscape_data/shrubs_27cellsr_sum.tif",
               "data/landscape_data/gramanoids_27cellsr_sum.tif", 
               "data/landscape_data/foodcrops_27cellsr_sum.tif", 
               "data/landscape_data/deciduous_27cellsr_sum.tif",
               "data/landscape_data/mixed_27cellsr_sum.tif",
               "data/landscape_data/evergreen_27cellsr_sum.tif", 
               "data/landscape_data/bottomland_27cellsr_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 2289

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

# Rename layers
names(layers) <- c("shrubs", "gramanoids", "foodcrops", "deciduous", "mixed", "evergreen", "bottomland", "water")


#### Extract covariates and used and available locations ####

## Reformat shapefiles
# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(layers))

# convert locations to spatvector
deer_v <- vect(deer, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)

# Extract distance values at used and available locations
dat_deer <- extract(layers, deer_v)

# Join extracted data back to location data frame
deer <- bind_cols(deer, dat_deer)

# correlation matrix
cor(deer[,c(7:14)])

# create key column
deer <- deer %>%
  unite("key", c("id", "burst"), sep="_", remove=FALSE)

# write file so we don't always have to wait for the rasters to do stuff
write_csv(deer, "output/deer_used_avail_covariates.csv")



