library(tidyverse)
library(terra)
library(sf)
library(exactextractr)

## Read location data
deer <- read_csv("data/location_data/deer_used_avail.csv")

# convert to sf object
deer <- st_as_sf(deer, coords=c("X", "Y"), crs=32616)

# Reproject to ACEA
deer <- st_transform(deer, crs=5070)

## Load rasters
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_buff <- vect("data/landscape_data/mississippi_ACEA_50kmbuffer.shp")

deer_buffer <- vect("data/landscape_data/deer_10km_buffers.shp")

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

library(lme4)
set.seed(1)

# Run model
m1 <- glmer(type ~ pctBottomland + pctDecid + pctEvergreen + pctWater + pctHerbaceous + pctCrops + (1|DeerID), data=deer, family=binomial(link = "logit"))

## Create empty raster
temp_rast <- rast(ext(ms_full), resolution=1000) # create template raster 1 km x 1 km

# Predict
pred_deer <- predict(layers, m1, type="response", re.form = NA)

# Mask to Mississippi
pred_deer <- mask(pred_deer, ms_full)

# Resample to 1 x 1 km for RAMAS
pred_deer <- resample(pred_deer, temp_rast)
pred_deer <- mask(pred_deer, ms_full)
pred_deer <- crop(pred_deer, ext(ms_full))
plot(pred_deer)

# Rescale to be between 0 and 1
pred_deer <- (pred_deer-minmax(pred_deer)[1])/(minmax(pred_deer)[2]-minmax(pred_deer)[1])
plot(pred_deer)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
pred_deer <- classify(pred_deer, m)

# writeRaster(pred_deer, "data/predictions/deer_glmm_rsf_FINALFINALFINALFINAL.tif", overwrite=TRUE)
# writeRaster(pred_deer, "data/predictions/deer_glmm_rsf_FINALFINALFINALFINAL.asc",NAflag=-9999, overwrite=TRUE)
# 

all_vals <- as.data.frame(pred_deer)
range(all_vals$lyr1)
mean(all_vals$lyr1)
quantile(all_vals$lyr1, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

## Extract RSF predictions at each location
deer_test <- deer %>% select(DeerID:Y)

dat_deer <- extract(pred_deer, deer_test[,c(4:5)])

deer_test$rsf <- dat_deer$lyr1

hist(deer_test$rsf[deer_test$type==1], xlim=c(0,1))
hist(deer_test$rsf[deer_test$type==0], xlim=c(0,1))

min(deer_test$rsf[deer_test$type==1], na.rm=TRUE)
mean(deer_test$rsf[deer_test$type==1], na.rm=TRUE)
median(deer_test$rsf[deer_test$type==1], na.rm=TRUE)

quantile(deer_test$rsf[deer_test$type==1], probs=c(0.25,0.5,0.75, 0.95))


thresh <- seq(0,1,by=0.01)
res <- matrix(NA, ncol=3, nrow=length(thresh))
res[,1] <- thresh
for (i in 1:length(thresh)) {
  
  actu <- factor(deer_test$type, levels=c(0,1))
  pred <- factor(if_else(deer_test$rsf >= thresh[i], 1, 0), levels=c(0,1))
  
  res[i,2] <- caret::sensitivity(pred, actu)
  res[i,3] <- caret::specificity(pred, actu)
}

res <- res %>% as.data.frame() 
names(res) <- c("Threshold", "Sensitivity", "Specificity")
res <- res %>% mutate(TSS = Sensitivity + Specificity - 1)

res %>% slice_max(TSS)


