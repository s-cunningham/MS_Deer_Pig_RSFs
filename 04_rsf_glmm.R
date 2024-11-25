library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(tidyterra)

options(scipen=999)
# source("00_functions.R")

#### Load rasters for prediction ####
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 

deer_buffer <- vect("data/landscape_data/deer_10km_buffers.shp")
pigs_buffer <- vect("data/landscape_data/PigsGPS_10kmBuffer.shp")

## Load rasters
rast_list <- c("data/landscape_data/distance_to_upland90m.tif",
               "data/landscape_data/distance_to_herbaceous90m.tif",
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

# Add to raster stack
layers <- c(layers, gpw, roads)

layers <- scale(layers)
names(layers) <- c("uplandforest", "grassland", "bottomlandhw", "water", "allforest", "ndvi_w", "ndvi_s", "foodcrops", "plantations", "streams", "tcc", "pop_density", "roads")

## Layers for each species model
deer_layers <- c(layers["ndvi_w"],  layers["foodcrops"], layers["streams"], layers["tcc"],  layers["grassland"]) #, layers["foodcrops"], layers["allforest"], 

pig_layers <- c(layers["bottomlandhw"], layers["uplandforest"], layers["foodcrops"], layers["roads"], layers["streams"])

#### Load data ####
## Load deer, add 5-fold label
deer <- read_csv("data/deer_usedavail_covariates.csv")
npts <- deer %>% filter(type==1) %>% group_by(DeerID) %>% count()
set.seed(11) # set seed to group roughly equal number of rows
npts$folds <- sample(1:5, length(unique(npts$DeerID)), replace=TRUE)
npts <- npts %>% select(-n)
deer <- left_join(deer, npts, by="DeerID") 

## Load pigs, add 5-fold label
pigs <- read_csv("data/pigs_usedavail_covariates.csv")
npts <- pigs %>% filter(type==1) %>% group_by(PigID) %>% count()
set.seed(8)
npts$folds <- sample(1:5, length(unique(npts$PigID)), replace=TRUE)
npts <- npts %>% select(-n)
pigs <- left_join(pigs, npts, by="PigID") 

set.seed(1)

### Deer model


## Cross-validation for deer
# Model 1
deer_layers <- c(layers["ndvi_w"],  layers["foodcrops"], layers["streams"], layers["tcc"],  layers["grassland"]) #, layers["foodcrops"], layers["allforest"], 
results <- data.frame(fold=1:5, buffer=NA, all=NA)
for (k in 1:5) {
  
  m1 <- glmer(type ~ ndvi_w + foodcrops + streams + tcc + grassland + (1|DeerID), data=deer[deer$folds!=k,], family=binomial(link = "logit"))
  
  pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
  
  # Mask raster to buffered area around points
  pred_buff <- mask(pred_deer, deer_buffer)
  
  # pred_deer <- mask(pred_deer, ms_full)
  pred_deer <- pred_deer / minmax(pred_deer)[2]
  
  # Get test points
  pts <- deer[deer$folds==k,4:5]
  
  boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE, plot=FALSE)
  results[k,2] <- boyce_buff$Boyce
  boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE, plot=FALSE)
  results[k,3] <- boyce_all$Boyce
}
saveRDS(results, "data/predictions/deer_m1.RDS")

## Model 2
deer_layers <- c(layers["allforest"], layers["foodcrops"], layers["streams"], layers["grassland"]) #, layers["foodcrops"], layers["allforest"], 
results <- data.frame(fold=1:5, buffer=NA, all=NA)
for (k in 1:5) {
  
  m1 <- glmer(type ~ allforest + foodcrops + streams + grassland + (1|DeerID), data=deer[deer$folds!=k,], family=binomial(link = "logit"))
  
  pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
  
  # Mask raster to buffered area around points
  pred_buff <- mask(pred_deer, deer_buffer)
  
  # pred_deer <- mask(pred_deer, ms_full)
  pred_deer <- pred_deer / minmax(pred_deer)[2]
  
  # Get test points
  pts <- deer[deer$folds==k,4:5]
  
  boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE, plot=FALSE)
  results[k,2] <- boyce_buff$Boyce
  boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE, plot=FALSE)
  results[k,3] <- boyce_all$Boyce
}
saveRDS(results, "data/predictions/deer_m2.RDS")

# Model 3
deer_layers <- c(layers["allforest"], layers["ndvi_w"], layers["foodcrops"], layers["streams"], layers["grassland"]) #, layers["foodcrops"], layers["allforest"], 
results <- data.frame(fold=1:5, buffer=NA, all=NA)
for (k in 1:5) {
  
  m1 <- glmer(type ~ allforest + ndvi_w + foodcrops + streams + grassland + (1|DeerID), data=deer[deer$folds!=k,], family=binomial(link = "logit"))
  
  pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
  
  # Mask raster to buffered area around points
  pred_buff <- mask(pred_deer, deer_buffer)
  
  # pred_deer <- mask(pred_deer, ms_full)
  pred_deer <- pred_deer / minmax(pred_deer)[2]
  
  # Get test points
  pts <- deer[deer$folds==k,4:5]
  
  boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE, plot=FALSE)
  results[k,2] <- boyce_buff$Boyce
  boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE, plot=FALSE)
  results[k,3] <- boyce_all$Boyce
}
saveRDS(results, "data/predictions/deer_m3.RDS")


# Model 4
deer_layers <- c(layers["allforest"], layers["ndvi_w"], layers["foodcrops"], layers["water"], layers["grassland"]) #, layers["foodcrops"], layers["allforest"], 
results <- data.frame(fold=1:5, buffer=NA, all=NA)
for (k in 1:5) {
  
  m1 <- glmer(type ~ allforest + ndvi_w + foodcrops + water + grassland + (1|DeerID), data=deer[deer$folds!=k,], family=binomial(link = "logit"))
  
  pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
  
  # Mask raster to buffered area around points
  pred_buff <- mask(pred_deer, deer_buffer)
  
  # pred_deer <- mask(pred_deer, ms_full)
  pred_deer <- pred_deer / minmax(pred_deer)[2]
  
  # Get test points
  pts <- deer[deer$folds==k,4:5]
  
  boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE, plot=FALSE)
  results[k,2] <- boyce_buff$Boyce
  boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE, plot=FALSE)
  results[k,3] <- boyce_all$Boyce
}
saveRDS(results, "data/predictions/deer_m4.RDS")


# Model 5
deer_layers <- c(layers["bottomlandhw"], layers["uplandforest"], layers["foodcrops"], layers["grassland"]) #, layers["foodcrops"], layers["allforest"], 
results <- data.frame(fold=1:5, buffer=NA, all=NA)
for (k in 1:5) {
  
  m1 <- glmer(type ~ bottomlandhw + uplandforest + foodcrops + grassland + (1|DeerID), data=deer[deer$folds!=k,], family=binomial(link = "logit"))
  
  pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
  
  # Mask raster to buffered area around points
  pred_buff <- mask(pred_deer, deer_buffer)
  
  # pred_deer <- mask(pred_deer, ms_full)
  pred_deer <- pred_deer / minmax(pred_deer)[2]
  
  # Get test points
  pts <- deer[deer$folds==k,4:5]
  
  boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE, plot=FALSE)
  results[k,2] <- boyce_buff$Boyce
  boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE, plot=FALSE)
  results[k,3] <- boyce_all$Boyce
}
saveRDS(results, "data/predictions/deer_m5.RDS")


# Model 6
deer_layers <- c(layers["bottomlandhw"], layers["plantations"], layers["foodcrops"], layers["grassland"]) #, layers["foodcrops"], layers["allforest"], 
results <- data.frame(fold=1:5, buffer=NA, all=NA)
for (k in 1:5) {
  
  m1 <- glmer(type ~ bottomlandhw + plantations + foodcrops + grassland + (1|DeerID), data=deer[deer$folds!=k,], family=binomial(link = "logit"))
  
  pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
  
  # Mask raster to buffered area around points
  pred_buff <- mask(pred_deer, deer_buffer)
  
  # pred_deer <- mask(pred_deer, ms_full)
  pred_deer <- pred_deer / minmax(pred_deer)[2]
  
  # Get test points
  pts <- deer[deer$folds==k,4:5]
  
  boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE, plot=FALSE)
  results[k,2] <- boyce_buff$Boyce
  boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE, plot=FALSE)
  results[k,3] <- boyce_all$Boyce
}
saveRDS(results, "data/predictions/deer_m6.RDS")


# Model 7
deer_layers <- c(layers["allforest"], layers["plantations"], layers["foodcrops"], layers["ndvi_w"], layers["water"]) #, layers["foodcrops"], layers["allforest"], 
results <- data.frame(fold=1:5, buffer=NA, all=NA)
for (k in 1:5) {
  
  m1 <- glmer(type ~ allforest + plantations + foodcrops + ndvi_w + water + (1|DeerID), data=deer[deer$folds!=k,], family=binomial(link = "logit"))
  
  pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
  
  # Mask raster to buffered area around points
  pred_buff <- mask(pred_deer, deer_buffer)
  
  # pred_deer <- mask(pred_deer, ms_full)
  pred_deer <- pred_deer / minmax(pred_deer)[2]
  
  # Get test points
  pts <- deer[deer$folds==k,4:5]
  
  boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE, plot=FALSE)
  results[k,2] <- boyce_buff$Boyce
  boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE, plot=FALSE)
  results[k,3] <- boyce_all$Boyce
}
saveRDS(results, "data/predictions/deer_m7.RDS")




pred_deer <- predict(deer_layers, m1, type="response", re.form = NA)
pred_deer <- mask(pred_deer, ms_full)
pred_deer <- pred_deer/ 0.2603978330 
plot(pred_deer)

# Mask raster to buffered area around points
pred_buff <- mask(pred_deer, deer_buffer)

pts <- deer %>% filter(type==1) %>% select(X, Y)

boyce_buff <- modEvA::Boyce(obs=pts, pred=pred_buff, rm.dup.points=TRUE)

boyce_all <- modEvA::Boyce(obs=pts, pred=pred_deer, rm.dup.points=TRUE)

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
