library(tidyverse)
library(terra)
library(tidyterra)

theme_set(theme_bw())

#### Read in scaled rasters ####
layers <- list.files(path="data/landscape_data/scaled_rasters/", pattern=".tif", full.names=TRUE)
layers <- layers[str_detect(layers, "deer")]
layers <- rast(layers)

# Read in MS shapefile (to drop islands)
ms <- vect("data/landscape_data/mississippi_ACEA.shp")
ms <- project(ms, layers)

# Read in permanent water mask
water <- vect("data/landscape_data/perm_water_grth500000m2.shp")
water <- project(water, layers)

#### Predict across MS counties ####
## Read in coefficients
# Betas
betas <- read_csv("output/deer_mixed_effects_betas.csv") 
betas[betas$covariate=="allhardwoods",1] <- "hardwoods"
betas[betas$covariate=="I(water_dist^2)",1] <- "water2"

beta_hat <- readRDS("output/deer_mixed_effects_beta.RDS")
names(beta_hat) <- c("hardwoods", "gramanoids", "foodcrops", "shrubs", "developed", "water_dist", "water2")
  
# Covariance
vcov_beta <- readRDS("output/deer_mixed_effects_vcov.RDS")

# drop intercept
vcov_beta <- vcov_beta[-1,-1]

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
  # rsf = hardwoods + evergreen + herbwetl + shrubs + gramanoids + developed + dist_Swater + I(dist_Swater^2) + dist_Pwater + I(dist_Pwater^2)
  pred <- exp(betas$beta[1]*cty_layers[["hardwoods"]] +
                betas$beta[2]*cty_layers[["gramanoids"]] +
                betas$beta[3]*cty_layers[["foodcrops"]] +
                betas$beta[4]*cty_layers[["shrubs"]] +
                betas$beta[5]*cty_layers[["developed"]] +
                betas$beta[6]*cty_layers[["water_dist"]] +
                betas$beta[7]*(cty_layers[["water_dist"]]^2))
  
  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_mean.tif")
  
  writeRaster(pred, filename, overwrite=TRUE)
  
  #### Plotting standard error ####
  # Create the squared covariate
  x_sq <- cty_layers[["water_dist"]]^2
  names(x_sq) <- "water2"  # name it to match the model formula
  
  # Add it to the covariate stack
  cty_layers <- c(cty_layers, x_sq)
  
  # Prediction matrix (X_pred)
  # Get the raster values as a matrix (rows = cells, cols = covariates)
  X_pred <- values(cty_layers, mat = TRUE)
  
  # rearrange column order
  X_pred <- X_pred[,c("hardwoods", "gramanoids", "foodcrops", "shrubs", "developed", "water_dist", "water2")]
  # Efficient matrix math
  SE <- sqrt(rowSums((X_pred %*% vcov_beta) * X_pred))
  
  # Just need 1 layer to form a template
  template <- cty_layers[[1]]
  
  # Convert to raster
  se_raster <- template # template raster
  terra::values(se_raster) <- SE
  names(se_raster) <- "SE"
  
  # create filename & save raster
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_SE.tif")
  writeRaster(se_raster, filename, overwrite=TRUE)
  
  ## Get 95% CI
  lin_pred <- betas$beta[1]*cty_layers[["hardwoods"]] +
    betas$beta[2]*cty_layers[["gramanoids"]] +
    betas$beta[3]*cty_layers[["foodcrops"]] +
    betas$beta[4]*cty_layers[["shrubs"]] +
    betas$beta[5]*cty_layers[["developed"]] +
    betas$beta[6]*cty_layers[["water_dist"]] +
    betas$beta[7]*(cty_layers[["water_dist"]]^2)
  
  # Compute upper and lower bounds
  lci_raster <- exp(lin_pred - 1.96 * se_raster)
  uci_raster <- exp(lin_pred + 1.96 * se_raster)
  
  # Lower
  names(lci_raster) <- "lower"
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_LCI.tif")
  writeRaster(lci_raster, filename, overwrite=TRUE)
  
  # Upper
  names(uci_raster) <- "upper"
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, "_UCI.tif")
  writeRaster(uci_raster, filename, overwrite=TRUE)
  
}

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
pred <- crop(pred, ms)

# Remove water
pred <- mask(pred, water, inverse=TRUE)

# plot
plot(pred)

m <- c(minmax(pred)[1], 0, 0,
       0, minmax(pred)[2], 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(pred, rclmat, include.lowest=TRUE)
plot(rc1)

# update layer and variable names
varnames(pred) <- "DeerResourceSelection"
names(pred) <- "RSF"


#### Check Boyce ####
library(modEvA)

# Calculate continuous Boyce index
deer <- read_csv("output/deer_used_avail_covariates.csv")

deerSV <- vect(deer, geom=c("X", "Y"), crs=crs(pred))
pt_preds <- extract(pred, deerSV)$RSF

obs <- deer$case

Boyce(obs=obs, pred=pt_preds)

## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_SE.tif$", full.names=TRUE)

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
# se_rast <- clamp(se_rast, upper=10, values=TRUE)

# plot SE raster (clamped at 10)
writeRaster(se_rast, "results/predictions/deer_rsf_se_30m.tif", overwrite=TRUE)


#### Resample to 90m ####
# Load template raster (cells must be *exactly* the same width/length for RAMAS)
temp_rast <- rast("data/landscape_data/CDL2023_90mACEA_mask.tif")

# Resample
pred90 <- resample(pred, temp_rast)
pred90 <- crop(pred90, pred)

## Calculate cutoff values for core and moderate areas
# convert raster to table, removing NAs
pred90_vals <- values(pred90, dataframe=TRUE)
pred90_vals <- pred90_vals |>
  as_tibble() |>
  filter(!is.na(RSF))

zero_val <- (0 - minmax(pred)[1])/(minmax(pred)[2]-minmax(pred)[1])

# Filter values above 0
above0 <- pred90_vals |>
  filter(RSF > 0)

# Calculate 50th percentile
# quantile(above0$RSF, probs=0.5)

# Where does 50th percentile fall on the linear stretch?
core <- (quantile(above0$RSF, probs=0.5) - minmax(pred)[1])/(minmax(pred)[2]-minmax(pred)[1])
core

marginal <- core / 2
marginal

#### Lower bound CI ####
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_LCI.tif$", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
lci_rast <- mosaic(rsrc)

# rename layer
names(lci_rast) <- "LowerCI"

# take ln of map
lci_rast <- lci_rast + 0.000001
lci_rast <- log(lci_rast)

# Reclassify missing data to 0
m <- rbind(c(NA, minmax(lci_rast)[1]))
lci_rast <- classify(lci_rast, m)

# Remove islands
lci_rast <- mask(lci_rast, ms)
lci_rast <- crop(lci_rast, ms)

# Remove water
lci_rast <- mask(lci_rast, water, inverse=TRUE)

plot(lci_rast)

#### Upper bound CI ####
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern="_UCI.tif$", full.names=TRUE)

## Read all files in as rasters
rlist <- lapply(files, rast)

## Convert to SpatRasterCollection
rsrc <- sprc(rlist)

## mosaic
uci_rast <- mosaic(rsrc)

# rename layer
names(uci_rast) <- "UpperCI"

# take ln of map
uci_rast <- uci_rast + 0.000001
uci_rast <- log(uci_rast)

# Reclassify missing data to 0
m <- rbind(c(NA, minmax(uci_rast)[1]))
uci_rast <- classify(uci_rast, m)

# Remove islands
uci_rast <- mask(uci_rast, ms)
uci_rast <- crop(uci_rast, ms)

# Remove water
uci_rast <- mask(uci_rast, water, inverse=TRUE)

plot(uci_rast)

#### Set up to write rasters for RAMAS ####

# Check min/max
minmax(pred)
minmax(lci_rast)
minmax(uci_rast)

# linear stretch (so that values >=0)
predls <- (pred - minmax(pred)[1])/(minmax(pred)[2]-minmax(pred)[1])

stretch_with_mean_range <- function(r, mean_rast) {
  (r - minmax(mean_rast)[1]) / (minmax(mean_rast)[2] - minmax(mean_rast)[1])
}

# Apply to CI rasters
lci_stretched <- stretch_with_mean_range(lci_rast, pred)
uci_stretched <- stretch_with_mean_range(uci_rast, pred)
# make sure its really between 0 and 1 (uci)...make all values >1 equal to 1
uci_stretched[uci_stretched > 1] <- 1

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
varnames(pred90) <- "DeerResourceSelection"
names(pred90) <- "RSF"

varnames(predls) <- "DeerResourceSelection"
names(predls) <- "RSF"

## Conf. Int

# Resample confidence interval rasters
lci90 <- resample(lci_stretched, temp_rast)
lci90 <- crop(lci90, lci_stretched)

uci90 <- resample(uci_stretched, temp_rast)
uci90 <- crop(uci90, uci_stretched)

## Write rasters

# Save ASCIIs for RAMAS
writeRaster(pred90, "results/predictions/deer_rsf_predicted.asc", NAflag=-9999, overwrite=TRUE)
writeRaster(lci90, "results/predictions/deer_rsf_lowerCI.asc", NAflag=-9999, overwrite=TRUE)
writeRaster(uci90, "results/predictions/deer_rsf_upperCI.asc", NAflag=-9999, overwrite=TRUE)


# additionally save mean predictions as .tif for plotting (both resampled and original 30x30m)
writeRaster(pred90, "results/predictions/deer_rsf_predicted_90m.tif", overwrite=TRUE)
writeRaster(pred, "results/predictions/deer_rsf_predicted_30m.tif", overwrite=TRUE)
writeRaster(predls, "results/predictions/deer_rsf_predicted_linStr_30m.tif", overwrite=TRUE)
# plot binned rsf values (log scale)
writeRaster(rc1, "results/predictions/deer_rsf_bins_30m.tif", overwrite=TRUE)

# Lower CI
writeRaster(lci_rast, "results/predictions/deer_rsf_LCI_30m.tif", overwrite=TRUE)

# Upper CI
writeRaster(uci_rast, "results/predictions/deer_rsf_UCI_30m.tif", overwrite=TRUE)



