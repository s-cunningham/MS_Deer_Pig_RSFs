library(tidyverse)
library(terra)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
deer <- read_csv("output/deer_used_avail_locations.csv") %>%
  # Add column for weight
  mutate(weight=if_else(case==1, 1, 5000))

# Rasters
rast_list <- c("data/landscape_data/shrublands_180m_sum.tif",
               "data/landscape_data/othercrops_180m_sum.tif",
               "data/landscape_data/gramanoids_180m_sum.tif", 
               "data/landscape_data/barren_180m_sum.tif", 
               "data/landscape_data/developed_180m_sum.tif",
               "data/landscape_data/palatable_crops_180m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 113

# Forest rasters (accounting for 50% mixed decid/evergreen)
forest <- c("data/landscape_data/evergreenMODIFIED_180m_sum.tif",
            "data/landscape_data/allhardwoodsMODIFIED_180m_sum.tif")
forest <- rast(forest)

# Reclassify missing data to 0
forest <- classify(forest, m)

# Convert to % 
forest <- forest / 226

layers <- c(layers, forest)

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

# Center and scale continuous rasters
layers <- scale(layers)

# Rename layers
names(layers) <- c("shrubs", "othercrops", "gramanoids", "barren", "developed", "foodcrops", "evergreen", "deciduous", "water")


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
cor(deer[,c(8:16)])

# create key column
deer <- deer %>%
  unite("key", c("id", "burst"), sep="_", remove=FALSE)

# write file so we don't always have to wait for the rasters to do stuff
write_csv(deer, "output/deer_used_avail_covariates.csv")
