# ----------------------------------------------------------------------------------------------------------------------------
# Name: 
# Author: Stephanie Cunningham
# Objective: Take the coefficients from the RSF and predict across Mississippi
# Input: Proportional cover rasters, coefficients from RSF models
# Output: Rasters (.tif and/or .asc) of statewide RSF predictions
# Required packages: tidyverse, terra, tidyterra
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(tidyterra)

#### Read in rasters ####
# Rasters
rast_list <- c("data/landscape_data/evergreen_180m_sum.tif",
               "data/landscape_data/deciduous_180m_sum.tif",
               "data/landscape_data/mixed_180m_sum.tif",
               "data/landscape_data/allhardwoods_180m_sum.tif",
               "data/landscape_data/shrublands_180m_sum.tif",
               "data/landscape_data/othercrops_180m_sum.tif",
               "data/landscape_data/gramanoids_180m_sum.tif", 
               "data/landscape_data/bottomlandHW_180m_sum.tif",
               "data/landscape_data/herbwetlands_180_sum.tif",
               "data/landscape_data/palatable_crops_180m_sum.tif",
               "data/landscape_data/developed_180m_sum.tif",
               "data/landscape_data/water_180m_sum.tif") 

layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 113

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)

layers <- c(layers, water)

# Rename layers
names(layers) <- c("evergreen", "deciduous", "mixed", "allhardwoods", "shrubs", "othercrops",
                   "gramanoids", "bottomland", "herbwetl", "foodcrops", "developed", "water_pct", "water_dist")

# remove some that aren't in the model
dontneed <- c("evergreen", "mixed", "othercrops", "herbwetl", "bottomland", "deciduous", "water_pct")
layers <- subset(layers, dontneed, negate=TRUE)

# Center and scale based on values in deer data
deer <- read_csv("output/deer_used_avail_covariates.csv")
# Center and scale covariates
deer_cs <- scale(deer[,7:19])
lnames <- names(layers)
for (i in 1:length(lnames)) {
  # Subtract mean and divide by standard deviation
  layers[[lnames[i]]] <- (layers[[lnames[i]]] - attr(deer_cs,"scaled:center")[[lnames[i]]]) / attr(deer_cs,"scaled:scale")[[lnames[i]]]
}

# Save pig center & scaled
# Define the output file names.  
output_files <- paste0("data/landscape_data/scaled_rasters/deer_scaled_", lnames, ".tif")

# Write each layer to a separate file
writeRaster(layers, output_files, overwrite=TRUE)


#### Predict across MS counties ####
# *Uses wayyyy less memory than doing the whole state at once*

## Read in coefficients
# Betas
betas <- read_csv("output/deer_glm_betas.csv") %>%
  # Pivot so that covariates are in a in a column, and beta coefficients in another column
  pivot_longer(1:ncol(.), names_to="covariate", values_to="beta")
# Standard error
se <- read_csv("output/deer_glm_se.csv") %>%
  # Pivot so that covariates are in a in a column, and beta coefficients in another column
  pivot_longer(1:ncol(.), names_to="covariate", values_to="se")

# Join into a single data frame
betas <- left_join(betas, se, by=c("covariate"))

# caluclate confidence intervales
betas <- betas %>%
  mutate(uci = beta + (1.96*se),
         lci = beta - (1.96*se))

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
  # allhardwoods + gramanoids + shrubs + foodcrops + developed + water
  pred <- exp(betas$beta[1]*cty_layers[["allhardwoods"]] +
                betas$beta[2]*cty_layers[["gramanoids"]] +
                betas$beta[3]*cty_layers[["shrubs"]] +
                betas$beta[4]*cty_layers[["foodcrops"]] +
                betas$beta[5]*cty_layers[["developed"]] +
                betas$beta[6]*cty_layers[["water_dist"]]+
                betas$beta[6]*(cty_layers[["water_dist"]]^2))
  
  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_mean.tif")

  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)
  
  # something not right with math
  ## Lower CI prediction
  pred <- exp(betas$lci[1]*cty_layers[["allhardwoods"]] +
                betas$lci[2]*cty_layers[["gramanoids"]] +
                betas$lci[3]*cty_layers[["shrubs"]] +
                betas$lci[4]*cty_layers[["foodcrops"]] +
                betas$lci[5]*cty_layers[["developed"]] +
                betas$lci[6]*cty_layers[["water_dist"]]+
                betas$lci[6]*(cty_layers[["water_dist"]]^2))

  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_LCI.tif")

  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)

  ## Upper CI prediction
  pred <- exp(betas$uci[1]*cty_layers[["allhardwoods"]] +
                betas$uci[2]*cty_layers[["gramanoids"]] +
                betas$uci[3]*cty_layers[["shrubs"]] +
                betas$uci[4]*cty_layers[["foodcrops"]] +
                betas$uci[5]*cty_layers[["developed"]] +
                betas$uci[6]*cty_layers[["water_dist"]]+
                betas$uci[6]*(cty_layers[["water_dist"]]^2))

  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_UCI.tif")

  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)
}

# Read in MS shapefile (to drop islands)
ms <- vect("data/landscape_data/mississippi_ACEA.shp")
ms <- project(ms, pred)

# Read in permanent water mask
water <- vect("data/landscape_data/perm_water_grth500000m2.shp")
water <- project(water, pred)

#### Combine county rasters into full state ####
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_mean.tif$", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
pred <- mosaic(rsrc)

# rename layer
names(pred) <- "RSF"

# take ln of map
pred <- pred + 0.000001
pred <- log(pred)

# Reclassify missing data to 0
m <- rbind(c(NA, minmax(pred)[1]))
pred <- classify(pred, m)

# Remove islands
pred <- mask(pred, ms)

# Remove water
pred <- mask(pred, water, inverse=TRUE)

# plot
plot(pred)

## Calculate quantiles
global(pred, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)

m <- c(-13.815511, -7.950021, 1,
       -7.950021, -5.639721, 2,
       -5.639721, -3.40708, 3,
       -3.40708, -2.085063, 4,
       -2.085063, -1.112339, 5,
       -1.112339, -0.3190405, 6,
       -0.3190405, 0.3061447, 7,
        0.3061447, 0.768716, 8,
        0.768716, 1.16346, 9,
        1.16346, 1.817576, 10)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(pred, rclmat, include.lowest=TRUE)
plot(rc1)

names(rc1) <- "bin"

# Set up NLCD classes and class values
rc_values <- 1:10
rc_class <- c("1", "2", "3", "4","5", "6", "7","8", "9", "10")

# Add class names and numbers to the raster
levels(rc1) <- list(data.frame(bin = rc_values,
                               suitability = rc_class))
# Plot
ggplot() +
  geom_spatraster(data=rc1) +
  scale_fill_grass_d(
    palette = "plasma",
    alpha = 1,
    direction = -1,
    na.translate = FALSE) +
  theme_void()

## Resample to 90x90
# Load template raster (cells must be *exactly* the same width/length for RAMAS)
temp_rast <- rast("data/landscape_data/CDL2023_90mACEA_mask.tif")

# Resample
pred90 <- resample(pred, temp_rast)
pred90 <- crop(pred90, pred)

# Check new quantiles
global(pred90, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)

# update layer and variable names
varnames(pred) <- "DeerResourceSelection"
names(pred) <- "RSF"

varnames(pred90) <- "DeerResourceSelection"
names(pred90) <- "RSF"

## Write rasters
# Save ASCII for RAMAS
writeRaster(pred90, "results/predictions/deer_rsf_predicted.asc", NAflag=-9999, overwrite=TRUE)
# additionally save mean predictions as .tif for plotting (both resampled and original 30x30m)
writeRaster(pred90, "results/predictions/deer_rsf_predicted_90m.tif", overwrite=TRUE)
writeRaster(pred, "results/predictions/deer_rsf_predicted_30m.tif", overwrite=TRUE)
# plot binned rsf values (log scale)
writeRaster(rc1, "results/predictions/deer_rsf_bins_30m.tif", overwrite=TRUE)

#### lower CI ####
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_LCI.tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
lci <- mosaic(rsrc)

# rename layer
names(lci) <- "RSF"

# take ln of map
lci <- lci + 0.000001
lci <- log(lci)

# Reclassify missing data to 0
m <- rbind(c(NA, minmax(lci)[1]))
lci <- classify(lci, m)

# Remove islands
lci <- mask(lci, ms)

# Remove water
lci <- mask(lci, water, inverse=TRUE)

# plot
plot(lci)

## Calculate quantiles
global(lci, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)





#### upper CI ####
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_UCI.tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
uci <- mosaic(rsrc)

# rename layer
names(uci) <- "RSF"

# take ln of map
uci <- uci + 0.000001
uci <- log(uci)

# Reclassify missing data to 0
m <- rbind(c(NA, minmax(uci)[1]))
uci <- classify(uci, m)

# Remove islands
uci <- mask(uci, ms)

# Remove water
uci <- mask(uci, water, inverse=TRUE)

# plot
plot(uci)

## Calculate quantiles
global(uci, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)


