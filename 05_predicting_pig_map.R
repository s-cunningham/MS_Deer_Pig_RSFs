library(tidyverse)
library(terra)
library(tidyterra)

theme_set(theme_bw())

#### Read Data ####
# pig data
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

# Center and scale covariates
pigs_cs <- scale(pigs[,7:17])

# Rasters
rast_list <- c("data/landscape_data/evergreen_210m_sum.tif",
               "data/landscape_data/deciduous_210m_sum.tif",
               "data/landscape_data/mixed_210m_sum.tif",
               "data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/othercrops_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/bottomlandHW_210m_sum.tif",
               "data/landscape_data/herbwetlands_210_sum.tif",
               "data/landscape_data/palatable_crops_210m_sum.tif",
               "data/landscape_data/developed_210m_sum.tif",
               "data/landscape_data/water_210m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# read water
# water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
# water <- resample(water, layers)
# ext(water) <- ext(layers)
# 
# layers <- c(layers, water)

# Rename layers
names(layers) <- c("evergreen", "deciduous", "mixed", "shrubs", "othercrops", "gramanoids", "bottomland", "herbwetl", "foodcrops", "developed", "water")

# remove some that aren't in the model
dontneed <- c("evergreen", "mixed", "othercrops", "herbwetl")
layers <- subset(layers, dontneed, negate=TRUE)

# Center and scale based on values in pig data
lnames <- names(layers)
for (i in 1:length(lnames)) {
  # Subtract mean and divide by standard deviation
  layers[[lnames[i]]] <- (layers[[lnames[i]]] - attr(pigs_cs,"scaled:center")[[lnames[i]]]) / attr(pigs_cs,"scaled:scale")[[lnames[i]]]
}

# Save pig center & scaled
# Define the output file names.  
output_files <- paste0("data/landscape_data/scaled_rasters/pig_scaled_", lnames, ".tif")

# Write each layer to a separate file
writeRaster(layers, output_files, overwrite=TRUE)


#### Predict across MS counties ####
## Read in coefficients
# Betas
betas <- read_csv("output/pigs_glm_betas.csv") %>%
  # Pivot so that covariates are in a in a column, and beta coefficients in another column
  pivot_longer(1:ncol(.), names_to="covariate", values_to="beta")
# Standard error
se <- read_csv("output/pigs_glm_se.csv") %>%
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

  # Predict (w(x) = exp(x*beta))
  # bottomland + deciduous + evergreen + shrubs + foodcrops + water + I(water^2)
  pred <- exp(betas$beta[1]*cty_layers[["bottomland"]] +
              betas$beta[2]*cty_layers[["shrubs"]] +
              betas$beta[3]*cty_layers[["developed"]] +
              betas$beta[4]*cty_layers[["water"]])
  
  # create filename
  filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, ".tif")
  
  writeRaster(pred, filename, overwrite=TRUE)
  
  ## Lower CI
  # pred <- exp(betas$lci[1]*cty_layers[["bottomland"]] +
  #               betas$lci[2]*cty_layers[["shrubs"]] +
  #               betas$lci[3]*cty_layers[["deciduous"]] +
  #               betas$lci[4]*cty_layers[["evergreen"]] +
  #               betas$lci[5]*cty_layers[["foodcrops"]] +
  #               betas$lci[6]*cty_layers[["water"]] +
  #               betas$lci[7]*(cty_layers[["water"]]^2))
  # 
  # # create filename
  # filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, "_LCI.tif")
  # 
  # # Export county prediction raster
  # writeRaster(pred, filename, overwrite=TRUE)
  # 
  # 
  # ## Upper CI
  # pred <- exp(betas$uci[1]*cty_layers[["bottomland"]] +
  #               betas$uci[2]*cty_layers[["shrubs"]] +
  #               betas$uci[3]*cty_layers[["deciduous"]] +
  #               betas$uci[4]*cty_layers[["evergreen"]] +
  #               betas$uci[5]*cty_layers[["foodcrops"]] +
  #               betas$uci[6]*cty_layers[["water"]] +
  #               betas$uci[7]*(cty_layers[["water"]]^2))
  # 
  # # create filename
  # filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, "_UCI.tif")
  # 
  # # Export county prediction raster
  # writeRaster(pred, filename, overwrite=TRUE)
  
}

#### Combine county rasters into full state ####
## List county rasters
files <- list.files(path="output/pig_county_preds/", pattern=".tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
pred <- mosaic(rsrc)

# rename layer
names(pred) <- "RSF"

# take ln of map
pred <- log(pred)

# plot
plot(pred)

## lower CI
## List county rasters
files <- list.files(path="output/pig_county_preds/", pattern="_LCI.tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
lci <- mosaic(rsrc)

# rename layer
names(lci) <- "RSF"

# take ln of map
lci <- log(lci)

# some values are -Inf...we don't want that, will replace with another low value
hist(lci)
# Reclassify missing data to 0
m <- rbind(c(-Inf, -103.27893))
lci <- classify(lci, m)

# plot
plot(lci)

# upper CI
## List county rasters
files <- list.files(path="output/pig_county_preds/", pattern="_UCI.tif", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
uci <- mosaic(rsrc)

# rename layer
names(uci) <- "RSF"

# take ln of map
uci <- log(uci)

# plot
plot(uci)

### 

# Check min/max of all three rasters
minmax(pred)
minmax(lci)
minmax(uci)

# Get quantile bins to use as thresholds
global(test_pred, fun=quantile,probs=seq(0.1, 1, 0.1), na.rm=TRUE)

global(lci, fun=quantile,probs=seq(0.1, 1, 0.1), na.rm=TRUE)

global(uci, fun=quantile,probs=seq(0.1, 1, 0.1), na.rm=TRUE)


## Need to resample for RAMAS (30 m cells going to have too many rows)
# Load template raster (cells must be *exactly* the same width/length for RAMAS)
temp_rast <- rast("data/landscape_data/CDL2023_90mACEA_mask.tif")

# Resample
pigs_pred <- resample(pred, temp_rast)
pigs_pred <- crop(pigs_pred, pred)

# Mask lakes
lakes <- vect("data/landscape_data/ne_ms_lakes.shp")
lakes <- project(lakes, pred)
pigs_pred <- mask(pigs_pred, lakes, inverse=TRUE)

# Reclassify missing data to 0
m <- rbind(c(NA, minmax(pred)[1]))
pigs_pred <- classify(pigs_pred, m)

plot(pigs_pred)



## Write rasters
writeRaster(pigs_pred, "results/predictions/pigs_rsf_predictedLn.asc", NAflag=-9999, overwrite=TRUE)
# additionally save mean predictions as .tif for plotting (both resampled and original 30x30m)
writeRaster(pigs_pred, "results/predictions/pigs_rsf_predictedLn_90m.tif", overwrite=TRUE)
writeRaster(pred, "results/predictions/pigs_rsf_predictedLn_30m.tif", overwrite=TRUE)

