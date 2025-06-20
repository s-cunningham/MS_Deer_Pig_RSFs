library(tidyverse)
library(terra)
library(tidyterra)


# get names of scaled rasters
files <- list.files("data/landscape_data/scaled_rasters/", pattern=".tif", full.names=TRUE)

# list layers
n <- c("deer_hardwoods", "deer_developed", "deer_crops", "deer_graminoids", "deer_shrubs", "deer_water", 
       "pigs_developed", "pigs_water", "pigs_graminoids", "pigs_hardwoods", "pigs_shrubs")

# read raster, rename layers
cov_rast <- rast(files)
names(cov_rast) <- n

## Clip to state outline
# Read in MS shapefile (to drop islands)
ms <- vect("data/landscape_data/mississippi_ACEA.shp")
ms <- project(ms, cov_rast)

# Read in permanent water mask
water <- vect("data/landscape_data/perm_water_grth500000m2.shp")
water <- project(water, cov_rast)

# Remove islands
cov_rast <- mask(cov_rast, ms)

# Remove water
cov_rast <- mask(cov_rast, water, inverse=TRUE)
names(cov_rast) <- n

plot(cov_rast[[6]])

# Replace all values >10 with 10
cov_rast <- clamp(cov_rast, lower=-5, upper=5, values = TRUE)
cov_rast <- crop(cov_rast, ms)


extrap_mask <- sum(abs(cov_stack) > 3, na.rm = TRUE) > 0