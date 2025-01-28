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
lakes <- vect("data/landscape_data/ne_ms_lakes.shp")
lakes <- project(lakes, ms_full)

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
# write_csv(deer, "data/deer_usedavail_HRcovariates.csv")

library(lme4)
set.seed(1)

## Run model
m1 <- glmer(type ~ pctBottomland + pctDecid + pctEvergreen + pctWater + pctHerbaceous + pctCrops + (1|DeerID), data=deer, family=binomial(link = "logit"))

## Create empty raster
temp_rast <- rast(ext(ms_full), resolution=1000) # create template raster 1 km x 1 km
temp_pts <- as.points(temp_rast, values=FALSE, na.rm=TRUE)
crs(temp_pts) <- crs(ms_full) # Define projection

# Predict
layers <- c(layers["pctBottomland"],  layers["pctDecid"], layers["pctEvergreen"], layers["pctWater"], layers["pctHerbaceous"], layers["pctCrops"])

#### Get SE of predictions
# Load points of raster cells
pts <- vect("data/landscape_data/ms_masked_90m_pts.shp")
pts <- project(pts, crs(layers))

# Extract raster values at points (i.e., put the raster values into a data frame)
pts_extr <- extract(layers, pts)
pts_extr <- pts_extr %>% select(-ID)

# Predict standard errors. Change type to "link" so that everything is on the same (log-odds) scale, 
# instead of generating predicted probabilities (we will do this conversion manually)
deer_sd <- predict(m1, newdata=pts_extr, type="link", se.fit=TRUE, re.form=NA)

# Create tibble for calculations
deer_sd <- data.frame(fit=deer_sd$fit, se.fit=deer_sd$se.fit)

# Calculate confidence intervals
deer_sd <- deer_sd %>% as_tibble() %>%
  mutate(lci=boot::inv.logit(fit - (1.96 * se.fit)),
         uci=boot::inv.logit(fit + (1.96 * se.fit))) %>%
  mutate(fit=boot::inv.logit(fit))

# Add values to pts SpatVector object. This takes up a lot of memory
pts$fit <- deer_sd$fit
pts$lci <- deer_sd$lci
pts$uci <- deer_sd$uci

# Write shapefile (slow...)
writeVector(pts, "data/landscape_data/deer_pts_se.shp", overwrite=TRUE)
# go open the point file in ArcGIS and convert it to a raster

## Convert points to raster. this is pretty memory-intensive (i.e., I hit 45 GB on my 64-GB machine)
deer_uci <- rasterize(pts, layers["pctBottomland"], field="uci", fun="mean")
deer_lci <- rasterize(pts, layers["pctBottomland"], field="lci", fun="mean")
deer_mean <- rasterize(pts, layers["pctBottomland"], field="fit", fun="mean")

## Read back in raster from points to do some final clipping/cleaning
# deer_uci <- rast("data/predictions/deer_glmm_rsf_UpperCI.tif")
# deer_lci <- rast("data/predictions/deer_glmm_rsf_LowerCI.tif")
# deer_mean <- rast("data/predictions/deer_glmm_rsf_mean.tif")

# Crop rows of NAs
deer_uci <- mask(deer_uci, ms_full)
deer_lci <- mask(deer_lci, ms_full)
deer_mean <- mask(deer_mean, ms_full)

## Rescale
# find min and max values
max_hsi <- max(minmax(deer_uci)[2], minmax(deer_lci)[2], minmax(deer_mean)[2])
min_hsi <- min(minmax(deer_uci)[1], minmax(deer_lci)[1], minmax(deer_mean)[1])

# apply scaling
deer_uci <- (deer_uci-min_hsi)/(max_hsi-min_hsi)
deer_lci <- (deer_lci-min_hsi)/(max_hsi-min_hsi)
deer_mean <- (deer_mean-min_hsi)/(max_hsi-min_hsi)

# Mask lakes
deer_uci <- mask(deer_uci, lakes, inverse=TRUE)
deer_lci <- mask(deer_lci, lakes, inverse=TRUE)
deer_mean <- mask(deer_mean, lakes, inverse=TRUE)

# Resample to 1x1 km
deer_uci <- resample(deer_uci, temp_rast)
deer_lci <- resample(deer_lci, temp_rast)
deer_mean <- resample(deer_mean, temp_rast)

# Get rid of extra space around MS
deer_uci <- crop(deer_uci, ext(ms_full))
deer_lci <- crop(deer_lci, ext(ms_full))
deer_mean <- crop(deer_mean, ext(ms_full))

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
deer_uci <- classify(deer_uci, m)
deer_lci <- classify(deer_lci, m)
deer_mean <- classify(deer_mean, m)

## Save LCI and UCI rasters as .asc
writeRaster(deer_uci, "data/predictions/deer_glmm_rsf_UpperCI.asc", NAflag=-9999, overwrite=TRUE)
writeRaster(deer_lci, "data/predictions/deer_glmm_rsf_LowerCI.asc", NAflag=-9999, overwrite=TRUE)
writeRaster(deer_mean, "data/predictions/deer_glmm_rsf_FINALFINALFINALFINAL.asc", NAflag=-9999, overwrite=TRUE)
# additionally save mean predictions as .tif for plotting
writeRaster(deer_mean, "data/predictions/deer_glmm_rsf_FINALFINALFINALFINAL.tif", overwrite=TRUE)

#### Check range of values and points ####
all_vals <- as.data.frame(deer_mean)
range(all_vals$mean)
mean(all_vals$mean)
quantile(all_vals$mean, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

## Extract RSF predictions at each location
deer_test <- deer %>% select(DeerID:Y)

dat_deer <- extract(deer_mean, deer_test[,c(4:5)])

deer_test$rsf <- dat_deer$mean

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


