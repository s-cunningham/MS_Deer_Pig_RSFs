library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(tidyterra)
library(ggspatial)

options(scipen=999)

#### Load rasters for prediction ####
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_buff <- vect("data/landscape_data/mississippi_ACEA_50kmbuffer.shp")

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

## foodcrops
food <- rast("data/landscape_data/distance_to_foodcrops90m.tif")
food <- project(food, crs(layers))
layers <- c(layers, food)

# USA streams rivers (ESRI)
streams <- rast("data/landscape_data/distance_to_rivers_streams.tif")
streams <- resample(streams, food, method="bilinear")
layers <- c(layers, streams)

gpw <- rast("data/landscape_data/gpw_PopDensity_90m50km.tif")
gpw <- project(gpw, crs(layers))
gpw <- resample(gpw, layers)
gpw <- mask(gpw, layers[[1]])
ext(gpw) <- ext(layers)
layers <- c(layers, gpw)

roads <- rast("data/landscape_data/distance_to_roads50kmbuff.tif")
roads <- project(roads, crs(layers))
ext(roads) <- ext(layers)
roads <- resample(roads, layers)
roads <- mask(roads, layers[[1]])
layers <- c(layers, roads)

layers <- scale(layers)
names(layers) <- c("uplandforest", "decidmixed", "grassland", "evergreen", "bottomlandhw", "water", "allforest", "foodcrops", "streams", "pop_density", "roads")

## Load pig locations, add 5-fold label
pigs <- read_csv("data/pigs_usedavail_covariates.csv")

#### Pig map ####
pigs_rsf <- glmer(type ~ bottomlandhw + uplandforest + foodcrops + roads + streams + (1|PigID), data=pigs, family=binomial(link = "logit"))
# pigs_rsf <- glmer(type ~ pctDecidMixed + pctCrops + pctBottomland + roads + streams + (1|PigID), data=pigs, family=binomial(link = "logit"))

# Predict across Mississippi
pig_layers <- c(layers["bottomlandhw"],  layers["uplandforest"], layers["foodcrops"], layers["roads"], layers["streams"])
pred_pigs <- predict(pig_layers, pigs_rsf, type="response", re.form = NA)

# Mask to state
pred_pigs <- mask(pred_pigs, ms_full)
plot(pred_pigs)

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

#### Prediction uncertainty ####
## Get SE of predictions
# Load points of raster cells
pts <- vect("data/landscape_data/ms_masked_90m_pts.shp")
pts <- project(pts, crs(layers))

# Extract raster values at points (i.e., put the raster values into a data frame)
pts_extr <- extract(pig_layers, pts)
pts_extr <- pts_extr %>% select(-ID)

# Predict standard errors. Change type to "link" so that everything is on the same (log-odds) scale, 
# instead of generating predicted probabilities (we will do this conversion manually)
pigs_sd <- predict(pigs_rsf, newdata=pts_extr, type="link", se.fit=TRUE, re.form=NA)

# Create tibble for calculations
pigs_sd <- data.frame(fit=pigs_sd$fit, se.fit=pigs_sd$se.fit)

# Calculate confidence intervals
pigs_sd <- pigs_sd %>% as_tibble() %>%
  mutate(lci=boot::inv.logit(fit - (1.96 * se.fit)),
         uci=boot::inv.logit(fit + (1.96 * se.fit))) %>%
  mutate(fit=boot::inv.logit(fit))

# Rescale confidence intervals based on rescaled predicted probabilities
pigs_sd <- pigs_sd %>% 
  mutate(fit=scales::rescale(fit, to=c(0,1)),
         lci=scales::rescale(lci, to=c(0,1)),
         uci=scales::rescale(uci, to=c(0,1))) 

pts$lci <- pigs_sd$lci
pts$uci <- pigs_sd$uci

# Write shapefile
writeVector(pts, "data/landscape_data/pigs_pts_se.shp", overwrite=TRUE)
# go open the point file in ArcGIS and convert it to a raster

## Read back in raster from points to do some final clipping/cleaning
pigs_uci <- rast("data/predictions/pigs_glmm_rsf_UpperCI.tif")
pigs_lci <- rast("data/predictions/pigs_glmm_rsf_LowerCI.tif")

# Resample and crop rows of NAs
pigs_uci <- mask(pigs_uci, ms_full)
pigs_lci <- mask(pigs_lci, ms_full)

pigs_uci <- resample(pigs_uci, temp_rast)
pigs_lci <- resample(pigs_lci, temp_rast)

pigs_uci <- crop(pigs_uci, ext(ms_full))
pigs_lci <- crop(pigs_lci, ext(ms_full))

plot(pigs_lci)
plot(pred_pigs)
plot(pigs_uci)

# Reclassify missing data to 0
pigs_uci <- classify(pigs_uci, m)
pigs_lci <- classify(pigs_lci, m)

## Save LCI and UCI rasters as .asc
writeRaster(pigs_uci, "data/predictions/pigs_glmm_rsf_UpperCI.asc",NAflag=-9999, overwrite=TRUE)
writeRaster(pigs_lci, "data/predictions/pigs_glmm_rsf_LowerCI.asc",NAflag=-9999, overwrite=TRUE)



#### Extract RSF predictions at each location ####
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







