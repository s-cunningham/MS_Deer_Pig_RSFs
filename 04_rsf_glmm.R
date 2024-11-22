library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(tidyterra)

#### Load rasters for prediction ####
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 

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

## Layers for each species model
deer_layers <- c(layers["allforest"], layers["plantations"], layers["foodcrops"], layers["streams"])

pig_layers <- c(layers["bottomlandhw"], layers["uplandforest"], layers["foodcrops"], layers["roads"], layers["streams"])

#### Load data ####
deer <- read_csv("data/deer_usedavail_covariates.csv")
pigs <- read_csv("data/pigs_usedavail_covariates.csv")

### Deer model

system.time(m1 <- glmer(type ~ tcc + plantations + foodcrops + streams + (1|DeerID), data=deer, family=binomial(link = "logit")))

pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)

pred_deer <- mask(pred_deer, ms_full)
pred_deer <- pred_deer / 0.35876536686 #0.95601763112

plot(pred_deer)

ggplot() +
  geom_spatraster(data=pred_deer) +
  scale_fill_whitebox_c(palette="viridi", na.value="transparent") +
  theme_void()
ggsave("figs/deer_glmm_rsf.svg")  


## Pig model
pigs_rsf <- glmer(type ~ bottomlandhw + uplandforest + foodcrops + roads + streams + (1|PigID), data=pigs, family=binomial(link = "logit"))

pred_pigs <- predict(pig_layers, pigs_rsf, type="response", re.form=NA)
pred_pigs <- mask(pred_pigs, ms_full)
pred_pigs <- pred_pigs/0.34445642786975 

ggplot() +
  geom_spatraster(data=pred_pigs) +
  scale_fill_whitebox_c(palette="viridi", na.value="transparent") +
  theme_void()
ggsave("figs/pigs_glmm_rsf.svg")  

pig_pts <- read_csv("data/location_data/pigs_used_avail.csv") %>% 
                filter(type==1) %>%
                st_as_sf(coords=c("X", "Y"), crs=32616) %>% 
                st_transform(crs=st_crs(pred_pigs))

pt_pred <- extract(pred_pigs, pig_pts)
hist(pt_pred$lyr1)


pigs4 <- read_csv("data/location_data/AllTrappingData.csv") %>%
  filter(Y>30 & Y<36 & no_hogs > 0) %>%
  rename(study=agency) %>% mutate(id=NA) %>%
  select(id, study, X, Y) %>% rename(Long=X, Lat=Y) %>%
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>%
  st_transform(crs=crs(pred_pigs))

trap <- extract(pred_pigs, pigs4)
hist(trap$lyr1)

## GBIF
pigs_gbif <- st_read("data/location_data/gbif_pigs.shp") %>% 
  st_transform(crs=4326) %>% st_coordinates() %>% as_tibble() %>% 
  rename(Long=X, Lat=Y) %>% mutate(study="gbif", id=NA) %>%
  select(id, study, Long, Lat) %>%
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>%
  st_transform(crs=crs(pred_pigs))

gbif <- extract(pred_pigs, pigs_gbif)
hist(gbif$lyr1)
