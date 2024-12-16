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
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_buff <- vect("data/landscape_data/mississippi_ACEA_50kmbuffer.shp")

deer_buffer <- vect("data/landscape_data/deer_10km_buffers.shp")
pigs_buffer <- vect("data/landscape_data/PigsGPS_10kmBuffer.shp")

## Load rasters
rast_list <- c("data/landscape_data/distance_to_upland90m.tif",
               "data/landscape_data/distance_to_decid_mixed90m.tif",
               "data/landscape_data/distance_to_herbaceous90m.tif",
               "data/landscape_data/distance_to_evergreen90m.tif",
               "data/landscape_data/distance_to_bottomlandhw90m.tif",
               "data/landscape_data/distance_to_openwater90m.tif",
               "data/landscape_data/distance_to_allforest90m.tif")
layers <- rast(rast_list)

## NDVI
# winter
ndvi_w <- rast("data/landscape_data/ndvi5year_jan1.tif")
ndvi_w <- project(ndvi_w, crs(layers))
ndvi_w  <- resample(ndvi_w, layers, method="bilinear")
layers <- c(layers, ndvi_w)

# summer
ndvi_s <- rast("data/landscape_data/ndvi5year_jun10.tif")
ndvi_s <- project(ndvi_s, crs(layers))
ndvi_s  <- resample(ndvi_s, layers, method="bilinear")
layers <- c(layers, ndvi_s)

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

## Soil properties
fc <- rast("data/landscape_data/fc_gNATSGO_MS50km_WGS84.tif")
fc <- project(fc, crs(layers))
fc <- resample(fc, layers, method="bilinear")
por <- rast("data/landscape_data/por_gNATSGO_MS50km_WGS84.tif") 
por <- project(por, crs(layers))
por <- resample(por, layers, method="bilinear")
awc <- rast("data/landscape_data/awc_gNATSGO_MS50km_WGS84.tif")
awc <- project(awc, crs(layers))
awc <- resample(awc, layers, method="bilinear")

layers <- c(layers, fc, por, awc)

## Net primary productivity
npp <- rast("data/landscape_data/2017-2023averageNPP_wgs84.tif")
npp <- project(npp, crs(layers))
npp <- resample(npp, layers, method="bilinear")

# Add to raster stack
layers <- c(layers, gpw, roads, npp)

## Percentage layers
pct_layers <- c("data/landscape_data/ms50km_pct_herbaceous_1pt5km.tif",
                "data/landscape_data/ms50km_pct_decidmixed_1pt5km.tif",
                "data/landscape_data/ms50km_pct_evergreen_1pt5km.tif",
                "data/landscape_data/ms50km_pct_bottomlandhw_1pt5km.tif",
                "data/landscape_data/ms50km_pct_medhighDevelopment_1pt5km.tif",
                "data/landscape_data/ms50km_pct_openWater_1pt5km.tif")
pct_layers <- rast(pct_layers)
pct_crops <- rast("data/landscape_data/ms50km_pct_foodcrops_1pt5km.tif")
pct_crops <- project(pct_crops, pct_layers)
pct_layers <- c(pct_layers, pct_crops)

pct_layers <- resample(pct_layers, layers, method="bilinear")

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
pct_layers <- classify(pct_layers, m)
pct_layers <- mask(pct_layers, ms_buff)


layers <- c(layers, pct_layers)


layers <- scale(layers)
names(layers) <- c("uplandforest", "decidmixed", "grassland", "evergreen", "bottomlandhw", "water", "allforest", "ndvi_w", "ndvi_s", "foodcrops",
                   "plantations", "streams", "tcc", "fc", "por", "awc", "pop_density", "roads", "npp", "pctHerbaceous", "pctDecidMixed", 
                   "pctEvergreen", "pctBottomland", "pctDevelopment", "pctWater", "pctCrops")


### high res pct rasters
# new_pct <- c("data/landscape_data/500m_pct_HerbShrubPast.tif",
#              "data/landscape_data/500m_pct_CultCrops.tif",
#              "data/landscape_data/500m_pct_Bottomlands.tif",
#              "data/landscape_data/500m_pct_Evergreen.tif",
#              "data/landscape_data/500m_pct_DecidMixed.tif",
#              "data/landscape_data/500m_pct_Water.tif")
# new_pct <- rast(new_pct)
# 
# m <- rbind(c(NA, 0))
# new_pct <- classify(new_pct, m)
# new_pct <- mask(new_pct, ms_buff)
# 
# new_pct <- scale(new_pct)
# names(new_pct) <- c("HerbShrubPast", "CultCrops", "Bottomlands", "Evergreen", "DecidMixed", "Water")
# 


# Extract distance values at used and available locations
dat_deer <- extract(layers, deer)
# dat_deer <- extract(new_pct, deer)
# dat_pigs <- extract(layers, pigs)

# Combine spatial data with points
deer_xy <- st_coordinates(deer)
deer <- st_drop_geometry(deer)
deer <- bind_cols(deer, dat_deer)
deer <- bind_cols(deer, deer_xy) %>% select(DeerID:ID,X,Y,uplandforest:pctCrops)
# deer <- bind_cols(deer, deer_xy) %>% select(DeerID:ID,X,Y,HerbShrubPast:Water)

pigs_xy <- st_coordinates(pigs)
pigs <- st_drop_geometry(pigs)
pigs <- bind_cols(pigs, dat_pigs)
pigs <- bind_cols(pigs, pigs_xy) %>% select(PigID:ID,X,Y,uplandforest:pctCrops)

# Check correlations
# cor(deer[,6:31])
cor(deer[,6:31])
cor(pigs[,6:31])

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

