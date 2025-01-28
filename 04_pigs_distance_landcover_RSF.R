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

lakes <- vect("data/landscape_data/ne_ms_lakes.shp")
lakes <- project(lakes, ms_full)

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

# create template 1 km x 1 km raster 
temp_rast <- rast(ext(ms_full), resolution=1000) 

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

pts$fit <- pigs_sd$fit
pts$lci <- pigs_sd$lci
pts$uci <- pigs_sd$uci
# pts$se <- pigs_sd$se.fit

# Write shapefile
writeVector(pts, "data/landscape_data/pigs_pts_se.shp", overwrite=TRUE)
# go open the point file in ArcGIS and convert it to a raster

## Convert points to raster. this is pretty memory-intensive
# The second argument is a raster that matches the resolution of the pts SpatVector (i.e., 90x90)
pigs_uci <- rasterize(pts, layers["bottomlandhw"], field="uci", fun="mean")
pigs_lci <- rasterize(pts, layers["bottomlandhw"], field="lci", fun="mean")
pigs_mean <- rasterize(pts, layers["bottomlandhw"], field="fit", fun="mean")

## Read back in raster from points to do some final clipping/cleaning
# pigs_uci <- rast("data/predictions/pigs_glmm_rsf_UpperCI.tif")
# pigs_lci <- rast("data/predictions/pigs_glmm_rsf_LowerCI.tif")
# pigs_mean <- rast("data/predictions/pigs_glmm_rsf_mean.tif")

# Crop rows of NAs
pigs_uci <- mask(pigs_uci, ms_full)
pigs_lci <- mask(pigs_lci, ms_full)
pigs_mean <- mask(pigs_mean, ms_full)

## Rescale
# find min and max values
max_hsi <- max(minmax(pigs_uci)[2], minmax(pigs_lci)[2], minmax(pigs_mean)[2])
min_hsi <- min(minmax(pigs_uci)[1], minmax(pigs_lci)[1], minmax(pigs_mean)[1])

# apply scaling
pigs_uci <- (pigs_uci-min_hsi)/(max_hsi-min_hsi)
pigs_lci <- (pigs_lci-min_hsi)/(max_hsi-min_hsi)
pigs_mean <- (pigs_mean-min_hsi)/(max_hsi-min_hsi)

# Mask lakes
pigs_uci <- mask(pigs_uci, lakes, inverse=TRUE)
pigs_lci <- mask(pigs_lci, lakes, inverse=TRUE)
pigs_mean <- mask(pigs_mean, lakes, inverse=TRUE)

# Resample to 1x1 km
pigs_uci <- resample(pigs_uci, temp_rast)
pigs_lci <- resample(pigs_lci, temp_rast)
pigs_mean <- resample(pigs_mean, temp_rast)

# Get rid of extra space around MS
pigs_uci <- crop(pigs_uci, ext(ms_full))
pigs_lci <- crop(pigs_lci, ext(ms_full))
pigs_mean <- crop(pigs_mean, ext(ms_full))

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
pigs_uci <- classify(pigs_uci, m)
pigs_lci <- classify(pigs_lci, m)
pigs_mean <- classify(pigs_mean, m)

## Save LCI and UCI rasters as .asc
writeRaster(pigs_uci, "data/predictions/pigs_glmm_rsf_UpperCI.asc",NAflag=-9999, overwrite=TRUE)
writeRaster(pigs_lci, "data/predictions/pigs_glmm_rsf_LowerCI.asc",NAflag=-9999, overwrite=TRUE)
writeRaster(pigs_mean, "data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.asc",NAflag=-9999, overwrite=TRUE)
# additionally save mean predictions as .tif for plotting
writeRaster(pigs_mean, "data/predictions/pigs_glmm_rsf_FINALFINALFINALFINAL.tif", overwrite=TRUE)

#### Extract RSF predictions at each location ####
pig_test <- pigs %>% select(PigID:Y)

dat_pig <- extract(pigs_mean, pig_test[,c(4:5)])

pig_test$rsf <- dat_pig$mean

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







