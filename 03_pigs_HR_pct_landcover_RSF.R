library(tidyverse)
library(terra)
library(sf)
library(lme4)

## Read location data
pigs <- read_csv("data/location_data/pigs_used_avail.csv")

# convert to sf object
pigs <- st_as_sf(pigs, coords=c("X", "Y"), crs=32616)

# Reproject to ACEA
pigs <- st_transform(pigs, crs=5070)

## Load rasters
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_buff <- vect("data/landscape_data/mississippi_ACEA_50kmbuffer.shp")

## Load rasters
rast_list <- c("data/landscape_data/14pt5km2_bottomlandSUM.tif",
               "data/landscape_data/14pt5km2_decidSUM.tif",
               "data/landscape_data/14pt5km2_evergreenSUM.tif",
               "data/landscape_data/14pt5km2_herbaceousSUM.tif",
               "data/landscape_data/14pt5km2_cropsSUM.tif")
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % of median home range
layers <- layers / 1793

names(layers) <- c("pctBottomland", "pctDecid", "pctEvergreen","pctHerbaceous","pctCrops")

# load distance to roads raster
roads <- rast("data/landscape_data/distance_to_roads50kmbuff.tif")
roads <- project(roads, crs(layers))
ext(roads) <- ext(layers)
roads <- resample(roads, layers)
roads <- mask(roads, layers[[1]])

# distance to streams
# USA streams rivers (ESRI)
streams <- rast("data/landscape_data/distance_to_rivers_streams.tif")
streams <- resample(streams, layers, method="bilinear")

# combine
layers <- c(layers, roads, streams)

names(layers)[6:7] <- c("roads", "streams")

# center & scale
layers <- scale(layers)

# Extract values from rasters
dat_pigs <- extract(layers, pigs)

# combine spatial data with points
pigs_xy <- st_coordinates(pigs)
pigs <- st_drop_geometry(pigs)
pigs <- bind_cols(pigs, dat_pigs)
pigs <- bind_cols(pigs, pigs_xy) %>% select(PigID:ID,X,Y,pctBottomland:streams)

# Check correlations
cor(pigs[,6:12])

# reorganize and create ID column
pigs <- pigs %>% separate(PigID, into=c("id", "study", "year"), sep="_") %>%
  select(-ID) %>%
  unite("PigID", 1:2, sep="_")

# Save file
write_csv(pigs, "data/pigs_usedavail_HRcovariates.csv")

#### Running RSF model
set.seed(1)

# Run model
m1 <- glmer(type ~ pctBottomland + pctDecid + pctEvergreen + pctCrops + roads + streams + (1|PigID), data=pigs, family=binomial(link = "logit"))

## Create empty raster
temp_rast <- rast(ext(ms_full), resolution=1000) # create template raster 1 km x 1 km

# Predict
pred_pigs <- predict(layers, m1, type="response", re.form = NA)

# Mask to Mississippi
pred_pigs <- mask(pred_pigs, ms_full)

# Resample to 1 x 1 km for RAMAS
pred_pigs <- resample(pred_pigs, temp_rast)
pred_pigs <- mask(pred_pigs, ms_full)
pred_pigs <- crop(pred_pigs, ext(ms_full))
plot(pred_pigs)

# Rescale to be between 0 and 1
pred_pigs <- (pred_pigs-minmax(pred_pigs)[1])/(minmax(pred_pigs)[2]-minmax(pred_pigs)[1])
plot(pred_pigs)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
pred_pigs <- classify(pred_pigs, m)

writeRaster(pred_pigs, "data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.tif", overwrite=TRUE)
writeRaster(pred_pigs, "data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.asc",NAflag=-9999, overwrite=TRUE)








