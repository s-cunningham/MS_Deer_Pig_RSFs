## Take the coefficients from the RSF and predict across Mississippi
## Because we are using 30 x 30 m cells, we have to go county by county

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
# Betas
betas <- read_csv("output/deer_glm_betas.csv") %>%
      # Pivot so that covariates are in a in a column, and beta coefficients in another column
      pivot_longer(1:7, names_to="covariate", values_to="beta")
# Standard error
se <- read_csv("output/deer_glm_se.csv") %>%
  # Pivot so that covariates are in a in a column, and beta coefficients in another column
  pivot_longer(1:7, names_to="covariate", values_to="se")

# Join into a single data frame
betas <- left_join(betas, se, by=c("covariate"))

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
  pred <- exp(betas$beta[1])*cty_layers[["deciduous"]] +
          exp(betas$beta[2])*cty_layers[["evergreen"]] +
          exp(betas$beta[3])*cty_layers[["gramanoids"]] +
          exp(betas$beta[4])*cty_layers[["shrubs"]] +
          exp(betas$beta[5])*cty_layers[["foodcrops"]] +
          exp(betas$beta[6])*cty_layers[["water"]] +
          exp(betas$beta[7])*(cty_layers[["water"]]^2)
  
  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_mean.tif")
  
  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)
  
  # something not right with math
  # # ## Lower CI prediction
  # pred <- exp(betas$se[1]*1.96-betas$beta[1])*cty_layers[["deciduous"]] +
  #         exp(betas$se[2]*1.96-betas$beta[2])*cty_layers[["evergreen"]] +
  #         exp(betas$se[3]*1.96-betas$beta[3])*cty_layers[["gramanoids"]] +
  #         exp(betas$se[4]*1.96-betas$beta[4])*cty_layers[["shrubs"]] +
  #         exp(betas$se[5]*1.96-betas$beta[5])*cty_layers[["foodcrops"]] +
  #         exp(betas$se[6]*1.96-betas$beta[6])*cty_layers[["water"]] +
  #         exp(betas$se[7]*1.96-betas$beta[7])*(cty_layers[["water"]]^2)
  # 
  # # create filename
  # filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_LCI.tif")
  # 
  # # Export county prediction raster
  # writeRaster(pred, filename, overwrite=TRUE)
  # 
  # # ## Upper CI prediction
  # pred <- exp(betas$se[1]*1.96+betas$beta[1])*cty_layers[["deciduous"]] +
  #         exp(betas$se[2]*1.96+betas$beta[2])*cty_layers[["evergreen"]] +
  #         exp(betas$se[3]*1.96+betas$beta[3])*cty_layers[["gramanoids"]] +
  #         exp(betas$se[4]*1.96+betas$beta[4])*cty_layers[["shrubs"]] +
  #         exp(betas$se[5]*1.96+betas$beta[5])*cty_layers[["foodcrops"]] +
  #         exp(betas$se[6]*1.96+betas$beta[6])*cty_layers[["water"]] +
  #         exp(betas$se[7]*1.96+betas$beta[7])*(cty_layers[["water"]]^2)
  # 
  # # create filename
  # filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_UCI.tif")
  # 
  # # Export county prediction raster
  # writeRaster(pred, filename, overwrite=TRUE)
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

# ## lower CI
# ## List county rasters
# files <- list.files(path="output/deer_county_preds/", pattern="_LCI.tif", full.names=TRUE)
# 
# ## Read all files in as rasters
# rlist <- lapply(files, rast)
# 
# ## Convert to SpatRasterCollection
# rsrc <- sprc(rlist)
# 
# ## mosaic
# lci <- mosaic(rsrc)
# 
# # rename layer
# names(lci) <- "RSF"
# 
# # plot
# plot(lci)
# 
# ## upper CI
# ## List county rasters
# files <- list.files(path="output/deer_county_preds/", pattern="_UCI.tif", full.names=TRUE)
# 
# ## Read all files in as rasters
# rlist <- lapply(files, rast)
# 
# ## Convert to SpatRasterCollection
# rsrc <- sprc(rlist)
# 
# ## mosaic
# uci <- mosaic(rsrc)
# 
# # rename layer
# names(uci) <- "RSF"
# 
# # plot
# plot(uci)

## Perform linear stretch
pred <- (pred - minmax(pred)[1])/(minmax(pred)[2] - minmax(pred)[1])
# pred <- (pred - minmax(pred)[1])/(minmax(pred)[2] - minmax(pred)[1])
# pred <- (pred - minmax(pred)[1])/(minmax(pred)[2] - minmax(pred)[1])

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
quantile(used$RSF, probs=c(0.05, 0.1,0.25, 0.5, 0.75,0.95, 1))


## Reclassify raster and export preliminary patches (just based on thresholds)
m <- c(0, 0.80, 0,
       0.80, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)

writeRaster(rc1, "results/predictions/deer_rsf_predicted_90m_core.tif", overwrite=TRUE)

m <- c(0, 0.40, 0,
       0.40, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)

writeRaster(rc1, "results/predictions/deer_rsf_predicted_90m_marginal.tif", overwrite=TRUE)


m <- c(0, 0.24, 0,
       0.24, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)

writeRaster(rc1, "results/predictions/deer_rsf_predicted_90m_highlymarginal.tif", overwrite=TRUE)
