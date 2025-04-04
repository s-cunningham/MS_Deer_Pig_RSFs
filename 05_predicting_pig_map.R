library(tidyverse)
library(terra)
library(recipes)

theme_set(theme_bw())

#### Custom Functions ####
# Function for scaling rasters
rast_scale <- function(raster, scaling_recipe) {
  
  # Extract the variable name from the raster
  var_name <- names(raster)
  
  # Extract the centering and scaling values from recipe
  center_value <- scaling_recipe$steps[[1]]$means[[var_name]]  # Assuming step_center was first
  scale_value <- scaling_recipe$steps[[2]]$sds[[var_name]]     # Assuming step_scale was second
  
  # Transform the raster
  r_centered <- raster - center_value
  r_scaled <- r_centered / scale_value
  
  return(r_scaled)
}


#### Read Data ####
# locations (used & available)
pigs <- read_csv("output/pigs_used_avail_locations.csv")

# Rasters
rast_list <- c("data/landscape_data/shrublands_210m_sum.tif",
               "data/landscape_data/othercrops_210m_sum.tif",
               "data/landscape_data/gramanoids_210m_sum.tif", 
               "data/landscape_data/bottomlandHW_210m_sum.tif",
               "data/landscape_data/decidmixed_210m_sum.tif",
               "data/landscape_data/evergreen_210m_sum.tif",
               "data/landscape_data/herbwetlands_210_sum.tif",
               "data/landscape_data/palatable_crops_210m_sum.tif") 
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
names(layers) <- c("shrubs", "othercrops", "gramanoids", "bottomland", "decidmixed", "evergreen", "herbwetlands", "foodcrops", "water")

#### Predict across MS counties ####
pigs <- read_csv("output/pigs_used_avail_covariates.csv") %>%
  select(case,key,shrubs:water)

## Center and scale covariates using the recipes package
# Determine the model formula and set which values to center + scale
scaling_recipe <- recipe(case ~ ., data=pigs) %>%
  step_center(all_numeric_predictors()) %>%
  step_scale(all_numeric_predictors()) %>%
  prep()

# Do the actual centering and scaling of the numeric predictors
pigs <- bake(scaling_recipe, pigs)

## Center/scale rasters using values from recipe
# Convert raster stack to list of rasters
layers <- as.list(layers)

# Run raster scaling function over list of rasters
layers <- lapply(layers, rast_scale, scaling_recipe=scaling_recipe)

# Convert centered/scaled raster list back to stack
layers <- rast(layers)

## Read in coefficients
betas <- read_csv("output/pigs_glm_betas.csv")

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
  # decidmixed + bottomland + foodcrops + herbwetlands + shrubs + gramanoids +  water + I(water^2)
  pred <- exp(betas$decidmixed[1]*cty_layers[["decidmixed"]] +
              betas$bottomland[1]*cty_layers[["bottomland"]] +
              betas$foodcrops[1]*cty_layers[["foodcrops"]] +
              betas$herbwetlands[1]*cty_layers[["herbwetlands"]] +
              betas$shrubs[1]*cty_layers[["shrubs"]] +
              betas$gramanoids[1]*cty_layers[["gramanoids"]] + 
              betas$water[1]*cty_layers[["water"]] +
              betas$water2[1]*(cty_layers[["water"]]^2))
  
  # create filename
  filename <- paste0("output/pig_county_preds/pig_pred_", split_counties[[i]]$COUNTYNAME, ".tif")
  
  writeRaster(pred, filename, overwrite=TRUE)
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

# plot
plot(pred)

## Perform linear stretch
pred <- (pred - minmax(pred)[1])/(minmax(pred)[2] - minmax(pred)[1])


#####

# ## from-to-becomes
# # classify the values into three groups 
# # all values >= 0 and <= 0.25 become 1, etc.
# m <- c(0, 0.7, 0,
#        0.7, 1, 1)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rc1 <- classify(pred, rclmat, include.lowest=TRUE)
# 
# plot(rc1)


#### Calculate threshold ####

## Extract RSF values at pigs locations
pigs <- read_csv("output/pigs_used_avail_covariates.csv") %>%
  # Drop covariate values
  select(X:case)

# convert locations to spatvector
pigs <- vect(pigs, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)         

# Extract RSF values at each pigs and available point
pigs_extr <- extract(pred, pigs)
pigs$RSF <- pigs_extr$RSF

# convert to data frame (will drop the X, Y coords, but thats ok)
pigs <- as.data.frame(pigs)

# ggplot(data=pigs) +
#   geom_histogram(aes(x=RSF)) +
#   facet_wrap(vars(case))
# 
# quantile(pigs$RSF[pigs$case==1], probs=c(0.25, 0.5, 0.75, 1))

## Calculate precision and recall (sensitivity) for case==1
# Create vector of thresholds
thresh <- seq(0.01,1,by=0.01)
res <- matrix(NA, ncol=3, nrow=length(thresh))
res[,1] <- thresh
for (i in 1:length(thresh)) {
  
  # Convert presence/available to factor
  actu <- factor(pigs$case, levels=c(0,1))
  # Convert RSF to binary based on threshold
  pred <- factor(if_else(pigs$RSF >= thresh[i], 1, 0), levels=c(0,1))
  
  cm <- caret::confusionMatrix(pred, actu, positive="1")
  
  # Calculate Sensitivity/precision
  res[i,2] <- unname(cm$byClass["Precision"])
  # Calculate recall
  res[i,3] <- cm$byClass["Recall"]
}

res <- res %>% as.data.frame() %>%
  rename(threshold=V1, Precision=V2, Recall=V3) %>%
  as_tibble() 



ggplot(res, aes(x=Recall, y=Precision)) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  geom_point() +
  theme_bw()

