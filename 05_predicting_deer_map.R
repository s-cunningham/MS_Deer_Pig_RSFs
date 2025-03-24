library(tidyverse)
library(terra)
library(tidyterra)

#### Read in rasters ####
# Rasters
rast_list <- c("data/landscape_data/shrublands_180m_sum.tif",
               "data/landscape_data/gramanoids_180m_sum.tif", 
               "data/landscape_data/palatable_crops_180m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 113

# Forest rasters (accounting for 50% mixed decid/evergreen)
forest <- c("data/landscape_data/evergreenMODIFIED_180m_sum.tif",
            "data/landscape_data/allhardwoodsMODIFIED_180m_sum.tif")
forest <- rast(forest)

# Reclassify missing data to 0
forest <- classify(forest, m)

# Convert to % 
forest <- forest / 226

layers <- c(layers, forest)

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

# Center and scale continuous rasters
layers <- scale(layers)

# Rename layers
names(layers) <- c("shrubs", "gramanoids", "foodcrops", "evergreen", "deciduous", "water")

#### Predict across MS counties ####
# *Uses wayyyy less memory than doing the whole state at once*

## Read in coefficients
betas <- read_csv("output/deer_lasso_betas_uncert.csv")

## Read county shapefile
counties <- vect("data/landscape_data/county_nrcs_a_ms.shp")
# Reproject to match rasters
counties <- project(counties, layers)

# Create a list where each element is a county
split_counties <- split(counties, "COUNTYNAME") 

# Loop over counties
for (i in 1:length(split_counties)) {
  
  # Subset just a single county
  cty <- split_counties[[i]]
  
  # Mask layers to county
  cty_layers <- crop(layers, cty, mask=TRUE)
  
  ## Mean prediction
  # Predict (w(x) = exp(x*beta))
  pred <- exp(betas$mean[1])*cty_layers[["deciduous"]] +
          exp(betas$mean[2])*cty_layers[["evergreen"]] +
          exp(betas$mean[3])*cty_layers[["gramanoids"]] +
          exp(betas$mean[4])*cty_layers[["shrubs"]] +
          exp(betas$mean[5])*cty_layers[["foodcrops"]] +
          exp(betas$mean[6])*cty_layers[["water"]] +
          exp(betas$mean[7])*(cty_layers[["water"]]^2)
  
  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_mean.tif")
  
  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)
  
  ## Lower CI prediction
  pred <- exp(betas$lci[1])*cty_layers[["deciduous"]] +
          exp(betas$lci[2])*cty_layers[["evergreen"]] +
          exp(betas$lci[3])*cty_layers[["gramanoids"]] +
          exp(betas$lci[4])*cty_layers[["shrubs"]] +
          exp(betas$lci[5])*cty_layers[["foodcrops"]] +
          exp(betas$lci[6])*cty_layers[["water"]] +
          exp(betas$lci[7])*(cty_layers[["water"]]^2)
  
  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_LCI.tif")
  
  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)
  
  ## Updder CI prediction
  pred <- exp(betas$uci[1])*cty_layers[["deciduous"]] +
          exp(betas$uci[2])*cty_layers[["evergreen"]] +
          exp(betas$uci[3])*cty_layers[["gramanoids"]] +
          exp(betas$uci[4])*cty_layers[["shrubs"]] +
          exp(betas$uci[5])*cty_layers[["foodcrops"]] +
          exp(betas$uci[6])*cty_layers[["water"]] +
          exp(betas$uci[7])*(cty_layers[["water"]]^2)
  
  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_UCI.tif")
  
  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)
}

#### Combine county rasters into full state ####
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_mean.tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
pred <- mosaic(rsrc)

# rename layer
names(pred) <- "RSF"

# plot
plot(pred)

## Perform linear stretch
pred <- (pred - minmax(pred)[1])/(minmax(pred)[2] - minmax(pred)[1])

#### Lower CI
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_LCI.tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
predLCI <- mosaic(rsrc)

# rename layer
names(predLCI) <- "RSF"

# plot
plot(predLCI)

#### Upper CI
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_UCI.tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
predUCI <- mosaic(rsrc)

# rename layer
names(predUCI) <- "RSF"

# plot
plot(predUCI)

minmax(pred)
minmax(predLCI)
minmax(predUCI)[2]

predLCI <- (predLCI - minmax(predLCI)[1])/(minmax(predLCI)[2] - minmax(predLCI)[1])
predUCI <- (predUCI - minmax(predUCI)[1])/(minmax(predUCI)[2] - minmax(predUCI)[1])




## Need to resample for RAMAS (30 m cells going to have too many rows)
# Load template raster (cells must be *exactly* the same width/length for RAMAS)
temp_rast <- rast("data/landscape_data/CDL2023_90mACEA_mask.tif")

# Resample
deer_pred <- resample(pred, temp_rast)
deer_pred <- crop(deer_pred, pred)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
deer_pred <- classify(deer_pred, m)

## Write rasters
writeRaster(deer_pred, "results/predictions/deer_rsf_predicted.asc", NAflag=-9999, overwrite=TRUE)
# additionally save mean predictions as .tif for plotting (both resampled and original 30x30m)
writeRaster(deer_pred, "results/predictions/deer_rsf_predicted_90m.tif", overwrite=TRUE)
writeRaster(pred, "results/predictions/deer_rsf_predicted_30m.tif", overwrite=TRUE)


## Classify out the core and marginal area from marginal and highly marginal, respectively
m <- c(0.82, 1, 0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)
# writeRaster(rc1, "results/predictions/deer_rsf_predicted90m_noCore.tif", overwrite=TRUE)
writeRaster(rc1, "results/predictions/deer_rsf_predicted90m_noCore.asc", NAflag=-9999, overwrite=TRUE)

m <- c(0.41, 1, 0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)
# writeRaster(rc1, "results/predictions/deer_rsf_predicted90m_noMarginalCore.tif", overwrite=TRUE)
writeRaster(rc1, "results/predictions/deer_rsf_predicted90m_noMarginalCore.asc", NAflag=-9999, overwrite=TRUE)

#### Calculate threshold ####

## Extract RSF values at deer locations
# Read in used/avail points
deer <- read_csv("output/deer_used_avail_covariates.csv") %>%
  # Drop covariate values
  select(X:case)

# Convert to SpatVector
deer_v <- vect(deer, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)
deer_v <- sf::st_as_sf(deer_v)

# subset to used locations
used <- deer %>% filter(case==1)

# convert locations to spatvector
used <- vect(used, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)         

# Extract RSF values at each deer and available point
used_extr <- extract(deer_pred, used)
used$RSF <- used_extr$RSF

# convert to data frame (will drop the X, Y coords, but thats ok)
used <- as.data.frame(used)
# drop extra columns and convert to tibble
used <- used %>% 
  select(key, case, RSF) %>% 
  as_tibble()

# Quantiles of used points
quantile(used$RSF, probs=c(0.1,0.25, 0.5, 0.75, 1))


## Reclassify raster and export preliminary patches (just based on thresholds)
m <- c(0, 0.82, 0,
       0.82, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)

writeRaster(rc1, "results/predictions/deer_rsf_predicted_90m_core.tif", overwrite=TRUE)

m <- c(0, 0.41, 0,
       0.41, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)

writeRaster(rc1, "results/predictions/deer_rsf_predicted_90m_marginal.tif", overwrite=TRUE)


m <- c(0, 0.246, 0,
       0.246, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)

writeRaster(rc1, "results/predictions/deer_rsf_predicted_90m_highlymarginal.tif", overwrite=TRUE)
