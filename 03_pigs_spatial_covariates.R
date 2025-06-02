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
  mutate(weight=if_else(case==1, 1, 5000)) %>%
  # drop 19212_Delta_4 - not HR, very exploratory
  filter(key!="19212_Delta_4")

# Plot the ratio of used to available points
pigs %>% group_by(key, case) %>% count() %>% 
  pivot_wider(names_from="case", values_from="n") %>%
  mutate(ratioUA=`1`/`0`) %>% 
  ggplot() + geom_density(aes(x=ratioUA))

# Read in rasters
rast_list <- c("data/landscape_data/evergreen_210m_sum.tif",
               "data/landscape_data/deciduous_210m_sum.tif",
               "data/landscape_data/mixed_210m_sum.tif",
               "data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/othercrops_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/bottomlandHW_210m_sum.tif",
               "data/landscape_data/herbwetlands_210_sum.tif",
               "data/landscape_data/palatable_crops_210m_sum.tif",
               "data/landscape_data/developed_210m_sum.tif",
               "data/landscape_data/allforestwoods_210m_sum.tif",
               "data/landscape_data/allhardwoods_210m_sum.tif",
               "data/landscape_data/water_210m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# Tree structure
# ch <- rast("data/landscape_data/LC23_CHavg_210m.tif")
# ch <- project(ch, layers)
# ext(ch) <- ext(ch)
# tcc <- rast("data/landscape_data/nlcd_tcc2021_ms50km_mean210m.tif")
# tcc <- project(tcc, layers)
# ext(tcc) <- ext(tcc)
# tcc_sd <- rast("data/landscape_data/nlcd_tcc2021_ms50km_sd210m.tif")
# tcc_sd <- project(tcc_sd, layers)
# ext(tcc_sd) <- ext(tcc_sd)

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)
seas_water <- rast("data/landscape_data/msGWD_seasonal_water_distance.tif")
seas_water <- project(seas_water, layers)
ext(seas_water) <- ext(seas_water)
perm_water <- rast("data/landscape_data/msGWD_permanent_water_distance.tif")
perm_water <- project(perm_water, layers)
ext(perm_water) <- ext(perm_water)

layers <- c(layers, water, seas_water, perm_water)

# Rename layers
names(layers) <- c("evergreen", "deciduous", "mixed", "shrubs", "othercrops", "gramanoids", "bottomland", "herbwetl", "foodcrops",
                   "developed", "allwoods", "hardwoods", "pct_water", "dist_water", "dist_Swater", "dist_Pwater")
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

cor(pigs[,7:22])

# write file so we don't always have to wait for the rasters to do stuff
write_csv(pigs, "output/pigs_used_avail_covariates.csv")



