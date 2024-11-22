library(tidyverse)
library(terra)
library(sf)
library(exactextractr)

## Read location data
deer <- read_csv("data/location_data/deer_used_avail.csv")
pigs <- read_csv("data/location_data/pigs_used_avail.csv")

# convert to sf object
deer <- st_as_sf(deer, coords=c("X", "Y"), crs=32616)
pigs <- st_as_sf(pigs, coords=c("X", "Y"), crs=32616)

# Reproject to ACEA
deer <- st_transform(deer, crs=5070)
pigs <- st_transform(pigs, crs=5070)

# create buffers
# deer_buff <- deer %>% st_buffer(100) %>% st_as_sf()

## Load rasters
rast_list <- c("data/landscape_data/distance_to_upland90m.tif",
               "data/landscape_data/distance_to_herbaceous90m.tif",
               "data/landscape_data/distance_to_bottomlandhw90m.tif",
               "data/landscape_data/distance_to_openwater90m.tif",
               "data/landscape_data/distance_to_allforest90m.tif")
layers <- rast(rast_list)

## foodcrops
food <- rast("data/landscape_data/distance_to_foodcrops90m.tif")
food <- project(food, crs(layers))
layers <- c(layers, food)

## plantations
plant <- rast("data/landscape_data/distance_to_plantation90m.tif")
plant <- project(plant, crs(layers))
layers <- c(layers, plant)

# USA streams rivers (ESRI)
streams <- rast("data/landscape_data/distance_to_rivers_streams.tif")
streams <- resample(streams, food, method="bilinear")
layers <- c(layers, streams)

# NLCD tree canopy cover
tcc <- rast("data/landscape_data/nlcd2021tcc_MWMean_90m.tif")  # EPSG 5070
tcc <- project(tcc, crs(layers))
layers <- c(layers, tcc)

# GPw, Roads, and streams need to be matched to extent and resampled to 90m
gpw <- rast("data/landscape_data/gpw_PopDensity_90m50km.tif")
gpw <- project(gpw, crs(layers))
gpw <- resample(gpw, layers)
gpw <- mask(gpw, layers[[1]])
ext(gpw) <- ext(layers)

roads <- rast("data/landscape_data/distance_to_roads50kmbuff.tif")
roads <- project(roads, crs(layers))
ext(roads) <- ext(layers)
roads <- resample(roads, layers)
roads <- mask(roads, layers[[1]])

# Add to raster stack
layers <- c(layers, gpw, roads)

layers <- scale(layers)
names(layers) <- c("uplandforest", "grassland", "bottomlandhw", "water", "allforest", "foodcrops", "plantations", "streams", "tcc", "pop_density", "roads")

# Extract distance values at used and available locations
dat_deer <- extract(layers, deer)
dat_pigs <- extract(layers, pigs)

# Combine spatial data with points
deer <- st_drop_geometry(deer)
deer <- bind_cols(deer, dat_deer)

pigs <- st_drop_geometry(pigs)
pigs <- bind_cols(pigs, dat_pigs)

# Check correlations
cor(deer[,4:13])
cor(pigs[,4:13])

## Deer
# Add a year column (split out ID)
deer <- deer %>% separate(DeerID, into=c("id", "study", "year"), sep="_") %>%
          select(-ID) %>%
          unite("DeerID", 1:2, sep="_")

# Save file
write_csv(deer, "data/deer_usedavail_covariates.csv")

## Pigs
pigs <- pigs %>% separate(PigID, into=c("id", "study", "year"), sep="_") %>%
  select(-ID) %>%
  unite("PigID", 1:2, sep="_")

# Save file
write_csv(pigs, "data/pigs_usedavail_covariates.csv")

