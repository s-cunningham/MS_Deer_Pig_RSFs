library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(tidyterra)
library(ggspatial)

options(scipen=999)
# source("00_functions.R")

#### Load rasters for prediction ####
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

## high res pct rasters
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


#### Load data ####


## Load pigs, add 5-fold label
pigs <- read_csv("data/pigs_usedavail_covariates.csv")
npts <- pigs %>% filter(type==1) %>% group_by(PigID) %>% count()
set.seed(8)
npts$folds <- sample(1:5, length(unique(npts$PigID)), replace=TRUE)
npts <- npts %>% select(-n)
pigs <- left_join(pigs, npts, by="PigID")

#### Pig map ####
pigs_rsf <- glmer(type ~ bottomlandhw + uplandforest + foodcrops + roads + streams + pop_density + (1|PigID), data=pigs, family=binomial(link = "logit"))
# pigs_rsf <- glmer(type ~ pctDecidMixed + pctCrops + pctBottomland + roads + streams + (1|PigID), data=pigs, family=binomial(link = "logit"))

# Predict across Mississippi
pig_layers <- c(layers["bottomlandhw"],  layers["uplandforest"], layers["foodcrops"], layers["roads"], layers["streams"], layers["pop_density"])
# pig_layers <- c(layers["pctDecidMixed"],  layers["pctCrops"], layers["pctBottomland"], layers["roads"], layers["streams"])
pred_pigs <- predict(pig_layers, pigs_rsf, type="response", re.form = NA)

# Mask to state
pred_pigs <- mask(pred_pigs, ms_full)

# create template 1 km x 1 km raster 
temp_rast <- rast(ext(ms_full), resolution=1000) 

# Resample and crop rows of NAs
pred_pigs <- resample(pred_pigs, temp_rast)
pred_pigs <- mask(pred_pigs, ms_full)
pred_pigs <- crop(pred_pigs, ext(ms_full))
plot(pred_pigs)

# Rescale to be between 0 and 1
# pred_pigs <- pred_pigs/minmax(pred_pigs)[2]
pred_pigs <- (pred_pigs-minmax(pred_pigs)[1])/(minmax(pred_pigs)[2]-minmax(pred_pigs)[1])
plot(pred_pigs)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
pred_pigs <- classify(pred_pigs, m)
plot(pred_pigs)

# Save tif for viewing in ArcGIS, and ASCII for reading into RAMAS
writeRaster(pred_pigs, "data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.tif", overwrite=TRUE)
writeRaster(pred_pigs, "data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.asc", NAflag=-9999, overwrite=TRUE)

## Extract RSF predictions at each location
pig_test <- pigs %>% select(PigID:Y)

dat_pig <- extract(pred_pigs, pig_test[,c(4:5)])

pig_test$rsf <- dat_pig$lyr1

plot(density(pig_test$rsf[pig_test$type==0]), xlim=c(0,1))
plot(density(pig_test$rsf[pig_test$type==1]), xlim=c(0,1))
hist(pig_test$rsf[pig_test$type==1], xlim=c(0,1))
hist(pig_test$rsf[pig_test$type==0], xlim=c(0,1))

min(pig_test$rsf[pig_test$type==1], na.rm=TRUE)
mean(pig_test$rsf[pig_test$type==1], na.rm=TRUE)

quantile(pig_test$rsf, probs=c(0.25, 0.5, 0.75))


thresh <- seq(0,1,by=0.01)
res <- matrix(NA, ncol=3, nrow=length(thresh))
res[,1] <- thresh
for (i in 1:length(thresh)) {
  
  actu <- factor(pig_test$type, levels=c(0,1))
  pred <- factor(if_else(pig_test$rsf >= thresh[i], 1, 0), levels=c(0,1))
  
  res[i,2] <- caret::sensitivity(pred, actu)
  res[i,3] <- caret::specificity(pred, actu)
}

res <- res %>% as.data.frame() 
names(res) <- c("Threshold", "Sensitivity", "Specificity")
res <- res %>% mutate(TSS = Sensitivity + Specificity - 1)

res %>% slice_max(TSS)



## Replace 0 with NA for plot
# m <- rbind(c(0,NA))
# pig_map <- classify(pred_pigs, m)
# deer_map <- classify(pred_deer, m)
pred_pigs <- mask(pred_pigs, ms_full)

pm <- ggplot() +
  geom_spatraster(data=pig_map) +
  scale_fill_whitebox_c(palette="viridi", na.value="transparent", limits=c(0,1), name="Habitat\nPreference") +
  ggtitle("Wild Pigs") +
  theme_bw() +
  theme(legend.position = "none")

dm <- ggplot() +
  geom_spatraster(data=deer_map) +
  scale_fill_whitebox_c(palette="viridi", na.value="transparent", limits=c(0,1), name="Habitat\nPreference") +
  ggtitle("White-tailed Deer") +
  theme_bw() +
  theme(legend.title=element_text(size=11),
        legend.text=element_text(size=10))

library(patchwork)
dm | pm
ggsave(filename="figs/rsf_plots.svg")











