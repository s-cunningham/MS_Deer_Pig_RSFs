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
pigs <- read_csv("output/pigs_used_avail_locations_indwt_30m.csv") %>%
  # Add column for weight
  # mutate(weight=if_else(case==1, 1, 5000)) %>%
  # drop 19212_Delta_4 - not HR, very exploratory
  filter(key!="19212_Delta_4") |>
  select(-avail, -used)


# Read in rasters
rast_list <- c("data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/palatable_crops_210m_sum.tif",
               "data/landscape_data/developed_210m_sum.tif",
               # "data/landscape_data/herbwetlands_210_sum.tif",
               # "data/landscape_data/bottomlandHW_210m_sum.tif",
               # "data/landscape_data/decidmixed_210m_sum.tif",
               "data/landscape_data/allhardwoods_210m_sum.tif") 
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

# Rename layers
names(layers) <- c("shrubs", "gramanoids", "foodcrops","developed", "hardwoods", "dist_water")

# Get mean and sd from entire covariate rasters
mean_vals <- global(layers, "mean", na.rm = TRUE)
sd_vals   <- global(layers, "sd", na.rm = TRUE)

# Apply to covariate rasters
cov_scaled <- (layers - mean_vals$mean) / sd_vals$sd

rast_cs <- bind_cols(mean_vals, sd_vals)
rast_cs <- rownames_to_column(rast_cs, var="layer")
write_csv(rast_cs, "output/pigs_raster_mean_sds.csv")

lnames <- names(layers)
output_files <- paste0("data/landscape_data/scaled_rasters/pigs_scaled_", lnames, ".tif")

# # Write each layer to a separate file
writeRaster(cov_scaled, output_files, overwrite=TRUE)

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

cor(pigs[,7:ncol(pigs)])

# write file so we don't always have to wait for the rasters to do stuff
write_csv(pigs, "output/pigs_used_avail_covariates_indwt_30m.csv")



