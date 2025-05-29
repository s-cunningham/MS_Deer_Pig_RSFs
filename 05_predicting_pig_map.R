library(tidyverse)
library(terra)
library(tidyterra)

theme_set(theme_bw())

#### Read Data ####
# pig data
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

# Center and scale covariates
pigs_cs <- scale(pigs[,7:22])

# Read in rasters
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
               "data/landscape_data/allforestwoods_210m_sum.tif",
               "data/landscape_data/allhardwoods_210m_sum.tif",
               "data/landscape_data/water_210m_sum.tif") 
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# read water
water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
water <- resample(water, layers)
ext(water) <- ext(layers)
seas_water <- rast("data/landscape_data/msGWD_seasonal_water_distance.tif")
seas_water <- project(seas_water, layers)
ext(seas_water) <- ext(seas_water)
perm_water <- rast("data/landscape_data/msGWD_permanent_water_distance.tif")
perm_water <- project(perm_water, layers)
ext(perm_water) <- ext(perm_water)

layers <- c(layers, water, seas_water, perm_water)

# Rename layers
names(layers) <- c("evergreen", "deciduous", "mixed", "shrubs", "othercrops", "gramanoids", "bottomland", "herbwetl", "foodcrops",
                   "developed", "allwoods", "hardwoods", "pct_water", "dist_water", "dist_Swater", "dist_Pwater")

# remove some that aren't in the model
dontneed <- c( "mixed", "deciduous", "allwoods", "foodcrops", "pct_water", "othercrops", "bottomland", "dist_Swater", "dist_Pwater", "pct_water")
layers <- subset(layers, dontneed, negate=TRUE)

# Center and scale based on values in pig data
lnames <- names(layers)
for (i in 1:length(lnames)) {
  # Subtract mean and divide by standard deviation
  layers[[lnames[i]]] <- (layers[[lnames[i]]] - attr(pigs_cs,"scaled:center")[[lnames[i]]]) / attr(pigs_cs,"scaled:scale")[[lnames[i]]]
}

# Save pig center & scaled
# Define the output file names.
# output_files <- paste0("data/landscape_data/scaled_rasters/pig_scaled_", lnames, ".tif")
# 
# # Write each layer to a separate file
# writeRaster(layers, output_files, overwrite=TRUE)


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
  # hardwoods + evergreen + herbwetl + shrubs + gramanoids + developed + dist_Swater + I(dist_Swater^2) + dist_Pwater + I(dist_Pwater^2)
  pred <- exp(betas$beta[1]*cty_layers[["hardwoods"]] +
              betas$beta[2]*cty_layers[["shrubs"]] +
              betas$beta[3]*cty_layers[["gramanoids"]] +
              betas$beta[4]*cty_layers[["developed"]] +
              betas$beta[5]*cty_layers[["dist_water"]] +
              betas$beta[6]*(cty_layers[["dist_water"]]^2))
  
  # create filename
  filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, "_mean.tif")
  
  writeRaster(pred, filename, overwrite=TRUE)
  
  ## Lower CI
  pred <- exp(betas$lci[1]*cty_layers[["hardwoods"]] +
                betas$lci[2]*cty_layers[["shrubs"]] +
                betas$lci[3]*cty_layers[["gramanoids"]] +
                betas$lci[4]*cty_layers[["developed"]] +
                betas$lci[5]*cty_layers[["dist_water"]] +
                betas$lci[6]*(cty_layers[["dist_water"]]^2))

  # create filename
  filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, "_LCI.tif")

  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)


  # ## Upper CI
  pred <- exp(betas$uci[1]*cty_layers[["hardwoods"]] +
                betas$uci[2]*cty_layers[["shrubs"]] +
                betas$uci[3]*cty_layers[["gramanoids"]] +
                betas$uci[4]*cty_layers[["developed"]] +
                betas$uci[5]*cty_layers[["dist_water"]] +
                betas$uci[6]*(cty_layers[["dist_water"]]^2))

  # create filename
  filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, "_UCI.tif")

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
files <- list.files(path="output/pig_county_preds/", pattern="_mean.tif$", full.names=TRUE)

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

m <- c(minmax(pred)[1], -6.802658, 1,
       -6.802658, -4.412018, 2,
       -4.412018, -3.009047, 3,
       -3.009047, -2.051441, 4,
       -2.051441, -1.376415, 5,
       -1.376415, -0.8999376, 6,
       -0.8999376, -0.4866923, 7,
       -0.4866923, 0.3623784, 8,
       0.3623784, 1.695516, 9,
       1.695516, minmax(pred)[2], 10)
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

ggplot() +
  geom_spatraster(data=rc1) +
  scale_fill_grass_d(
    palette = "viridis",
    alpha = 1,
    direction = 1,
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
writeRaster(pred90, "results/predictions/pigs_rsf_predicted.asc", NAflag=-9999, overwrite=TRUE)
# additionally save mean predictions as .tif for plotting (both resampled and original 30x30m)
writeRaster(pred90, "results/predictions/pigs_rsf_predicted_90m.tif", overwrite=TRUE)
writeRaster(pred, "results/predictions/pigs_rsf_predicted_30m.tif", overwrite=TRUE)
# plot binned rsf values (log scale)
writeRaster(rc1, "results/predictions/pigs_rsf_bins_30m.tif", overwrite=TRUE)
