# ----------------------------------------------------------------------------------------------------------------------------
# Name: Spatial Covariate Extraction for Wild Pigs
# Author: Stephanie Cunningham
# Objective: Extract spatial covariates from rasters (that had been generated in ArcGIS Pro) and add to used/avail locations
# Input: Summed cover rasters (from Cropland Data Layer), used/avail points from home range analysis
# Output: Data frame with covariate values appended to used/avail locations
# Required packages: tidyverse, terra
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(terra)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
pigs <- read_csv("output/pigs_used_avail_locations.csv") %>%
  # Add column for weight
  mutate(weight=if_else(case==1, 1, 5000))

# Rasters
rast_list <- c("data/landscape_data/evergreen_210m_sum.tif",
               "data/landscape_data/deciduous_210m_sum.tif",
               "data/landscape_data/mixed_210m_sum.tif",
               "data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/othercrops_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/bottomlandHW_210m_sum.tif",
               "data/landscape_data/herbwetlands_210_sum.tif",
               "data/landscape_data/palatable_crops_210m_sum.tif",
               "data/landscape_data/developed_210m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

# Center and scale continuous rasters
# layers <- scale(layers)

# Rename layers
names(layers) <- c("evergreen", "deciduous", "mixed", "shrubs", "othercrops", "gramanoids", "bottomland", "herbwetl", "foodcrops", "developed", "water")
# global(layers[["shrubs"]], fun="mean")
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

# Combine a couple classes
pigs <- pigs %>%
  mutate(short_veg=gramanoids+shrubs,
         allcrops=othercrops+foodcrops)

# Correlation matrix
ids <- unique(pigs$id)
quick.cor <- list()
for (i in 1:length(ids)) {
  temp <- pigs %>% filter(id==ids[i])
  
  c.mat <- cor(temp[,c(7:19)])
  
  high.cor <- matrix(0, nrow=13, ncol=13)
  high.cor[which(abs(c.mat)>0.7)] <- 1
  high.cor <- as.data.frame(high.cor)
  names(high.cor) <- names(temp[,7:19])
  high.cor$var <- names(temp[,7:19])
  
  quick.cor[[i]] <- high.cor
}

cor(pigs[,c(7:19)])

# create key column
# pigs <- pigs %>%
#   unite("key", c("id", "burst"), sep="_", remove=FALSE)




# write file so we don't always have to wait for the rasters to do stuff
write_csv(pigs, "output/pigs_used_avail_covariates.csv")
