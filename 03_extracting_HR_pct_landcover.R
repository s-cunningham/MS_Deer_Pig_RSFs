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
rast_list <- c("data/landscape_data/2pt3km2_bottomlandSUM.tif",
               "data/landscape_data/2pt3km2_decidSUM.tif",
               "data/landscape_data/2pt3km2_evergreenSUM.tif",
               "data/landscape_data/2pt3km2_developedSUM.tif",
               "data/landscape_data/2pt3km2_waterSUM.tif",
               "data/landscape_data/2pt3km2_herbaceousSUM.tif",
               "data/landscape_data/2pt3km2_cropsSUM.tif")
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % of median home range
layers <- layers / 317

layers <- scale(layers)
names(layers) <- c("pctBottomland", "pctDecid", "pctEvergreen","pctDeveloped","pctWater","pctHerbaceous","pctCrops")

# Extract distance values at used and available locations
dat_deer <- extract(layers, deer)

# Combine spatial data with points
deer_xy <- st_coordinates(deer)
deer <- st_drop_geometry(deer)
deer <- bind_cols(deer, dat_deer)
deer <- bind_cols(deer, deer_xy) %>% select(DeerID:ID,X,Y,pctBottomland:pctCrops)

cor(deer[,6:12])

## Deer
# Add a year column (split out ID)
deer <- deer %>% separate(DeerID, into=c("id", "study", "year"), sep="_") %>%
  select(-ID) %>%
  unite("DeerID", 1:2, sep="_")

# Save file
write_csv(deer, "data/deer_usedavail_HRcovariates.csv")



