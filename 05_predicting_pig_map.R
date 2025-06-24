library(tidyverse)
library(terra)
library(tidyterra)

theme_set(theme_bw())

#### Read Data ####
# pig data
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

# Read in rasters
rast_list <- c("data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/developed_210m_sum.tif",
               "data/landscape_data/allhardwoods_210m_sum.tif") 
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

layers <- c(layers, water)

# Rename layers
names(layers) <- c("shrubs", "gramanoids", "developed", "hardwoods", "dist_water")

# Read in MS shapefile (to drop islands)
ms <- vect("data/landscape_data/mississippi_ACEA.shp")
ms <- project(ms, layers)

# Remove islands
layers <- mask(layers, ms)

# Read in permanent water mask
water <- vect("data/landscape_data/perm_water_grth500000m2.shp")
water <- project(water, layers)

# Remove water
layers <- mask(layers, water, inverse=TRUE)
# crop to map extent
layers <- crop(layers, ms)

# put layer names back
names(layers) <- c("shrubs", "gramanoids", "developed", "hardwoods", "dist_water")

# Get mean and sd from entire covariate rasters
mean_vals <- global(layers, "mean", na.rm = TRUE)
sd_vals   <- global(layers, "sd", na.rm = TRUE)

# Apply to covariate rasters
layers <- (layers - mean_vals$mean) / sd_vals$sd

# Save deer center & scaled
# Define the output file names. 
lnames <- c("shrubs", "gramanoids", "developed", "hardwoods", "dist_water")
output_files <- paste0("data/landscape_data/scaled_rasters/pigs_scaled_", lnames, ".tif")

# # Write each layer to a separate file
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

# Extract coefficients and their variance-covariance matrix for SE plotting
coefs <- readRDS("output/pigs_vcov.rds")
beta_hat <- coefs[[1]]
vcov_mat <- coefs[[2]]

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
  
  #### Plotting standard error ####
  # Create the squared covariate
  x_sq <- cty_layers[["dist_water"]]^2
  names(x_sq) <- "water2"  # name it to match the model formula

  # Add it to the covariate stack
  cty_layers <- c(cty_layers, x_sq)
  
  # Prediction matrix (X_pred)
  # Get the raster values as a matrix (rows = cells, cols = covariates)
  X_pred <- values(cty_layers, mat = TRUE)
  
  eta_var <- rowSums((X_pred %*% vcov_mat * X_pred))  # efficient shortcut
  eta_se <- sqrt(eta_var)
  
  # Just need 1 layer to form a template
  template <- cty_layers[[1]]
  
  # Convert to raster
  r_se <- template; values(r_se) <- eta_se
  names(r_se) <- "SE"
  
  # create filename
  filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, "_SE.tif")
  
  writeRaster(r_se, filename, overwrite=TRUE)
}

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

m <- c(minmax(pred)[1], 0, 0,
       0, minmax(pred)[2], 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(pred, rclmat, include.lowest=TRUE)
plot(rc1)
  
# linear stretch (so that values >=0)
predls <- (pred - minmax(pred)[1])/(minmax(pred)[2]-minmax(pred)[1])

## Calculate quantiles
qrtls <- global(predls, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)

m <- c(0, qrtls[1,1], 1,
       qrtls[1,1], qrtls[1,2], 2,
       qrtls[1,2], qrtls[1,3], 3,
       qrtls[1,3], qrtls[1,4], 4,
       qrtls[1,4], qrtls[1,5], 5,
       qrtls[1,5], qrtls[1,6], 6,
       qrtls[1,6], qrtls[1,7], 7,
       qrtls[1,7], qrtls[1,8], 8,
       qrtls[1,8], qrtls[1,9], 9,
       qrtls[1,9], 1, 10)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(predls, rclmat, include.lowest=TRUE)
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
pred90 <- resample(predls, temp_rast)
pred90 <- crop(pred90, predls)

# Check new quantiles
global(pred90, quantile, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), na.rm=TRUE)

# update layer and variable names
varnames(pred) <- "PigResourceSelection"
names(pred) <- "RSF"

varnames(pred90) <- "PigResourceSelection"
names(pred90) <- "RSF"

varnames(predls) <- "PigResourceSelection"
names(predls) <- "RSF"

## Write rasters
# Save ASCII for RAMAS
writeRaster(pred90, "results/predictions/pigs_rsf_predicted.asc", NAflag=-9999, overwrite=TRUE)
# additionally save mean predictions as .tif for plotting (both resampled and original 30x30m)
writeRaster(pred90, "results/predictions/pigs_rsf_predicted_90m.tif", overwrite=TRUE)
writeRaster(pred, "results/predictions/pigs_rsf_predicted_30m.tif", overwrite=TRUE)
writeRaster(predls, "results/predictions/pigs_rsf_predicted_linStr_30m.tif", overwrite=TRUE)
# plot binned rsf values (log scale)
writeRaster(rc1, "results/predictions/pigs_rsf_bins_30m.tif", overwrite=TRUE)


library(modEvA)

# Calculate continuous Boyce index
pigs <- read_csv("output/pigs_used_avail_covariates.csv")

pigsSV <- vect(pigs, geom=c("X", "Y"), crs=crs(pred))
pt_preds <- extract(pred, pigsSV)$RSF

obs <- pigs$case

Boyce(obs=obs, pred=pt_preds)


#### Uncertainty ####

## List county rasters
files <- list.files(path="output/pig_county_preds/", pattern="_SE.tif$", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
se_rast <- mosaic(rsrc)

# rename layer
names(se_rast) <- "SE"

# Read in MS shapefile (to drop islands)
ms <- project(ms, se_rast)

# Remove islands
se_rast <- mask(se_rast, ms)

# Read in permanent water mask
water <- project(water, se_rast)

# Remove water
se_rast <- mask(se_rast, water, inverse=TRUE)

plot(se_rast)

# Replace all values >10 with 10
se_rast <- clamp(se_rast, upper=10, values=TRUE)

# plot SE raster (clamped at 10)
writeRaster(se_rast, "results/predictions/pigs_rsf_se_30m.tif", overwrite=TRUE)



