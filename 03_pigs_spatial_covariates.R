library(tidyverse)
library(terra)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
pigs <- read_csv("output/pigs_used_avail_locations.csv") %>%
            # Add column for weight
            mutate(weight=if_else(case==1, 1, 5000))

# Rasters
rast_list <- c("data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/othercrops_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/bottomlandHW_210m_sum.tif",
               # "data/landscape_data/herbwetlands_210_sum.tif",
               # "data/landscape_data/allwetlands_210_sum.tif",
               "data/landscape_data/palatable_crops_210m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# Forest rasters (accounting for 50% mixed decid/evergreen)
forest <- c("data/landscape_data/evergreenMODIFIED_210m_sum.tif",
            "data/landscape_data/allhardwoodsMODIFIED_210m_sum.tif")
forest <- rast(forest)

# Reclassify missing data to 0
forest <- classify(forest, m)

# Convert to % 
forest <- forest / 298

layers <- c(layers, forest)

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

# Center and scale continuous rasters
layers <- scale(layers)

# Rename layers
names(layers) <- c("shrubs", "othercrops", "gramanoids", "bottomland", "foodcrops", "evergreen", "deciduous", "water")
global(layers[["shrubs"]], fun="mean")
# global(layers[["shrubs"]], fun="sd")

#### Extract covariates and used and available locations ####

## reformat shapefiles
# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(layers))

# convert locations to spatvector
pigs_v <- vect(pigs, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)

# Extract distance values at used and available locations
dat_pigs <- extract(layers, pigs_v)

# Join extracted data back to location data frame
pigs <- bind_cols(pigs, dat_pigs)

# Correlation matrix
cor(pigs[,c(8:15)])

# create key column
pigs <- pigs %>%
  unite("key", c("id", "burst"), sep="_", remove=FALSE)

# write file so we don't always have to wait for the rasters to do stuff
write_csv(pigs, "output/pigs_used_avail_covariates.csv")
