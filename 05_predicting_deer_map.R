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
betas <- read_csv("output/deer_lasso_betas.csv")

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
  pred <- exp(betas$deciduous[1])*cty_layers[["deciduous"]] +
          exp(betas$evergreen[1])*cty_layers[["evergreen"]] +
          exp(betas$gramanoids[1])*cty_layers[["gramanoids"]] +
          exp(betas$shrubs[1])*cty_layers[["shrubs"]] +
          exp(betas$foodcrops[1])*cty_layers[["foodcrops"]] +
          exp(betas$water[1])*cty_layers[["water"]] +
          exp(betas$water2[1])*(cty_layers[["water"]]^2)
  
  # create filename
  filename <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME, ".tif")
  
  # Export county prediction raster
  writeRaster(pred, filename, overwrite=TRUE)
}

#### Combine county rasters into full state ####
## List county rasters
files <- list.files(path="output/deer_county_preds/", pattern=".tif", full.names=TRUE)

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

#### Calculate threshold ####
## Extract RSF values at deer locations
# Read in used/avail points
deer <- read_csv("output/deer_used_avail_covariates.csv") %>%
          # Drop covariate values
          select(X:case)

# Convert to SpatVector
deer_v <- vect(deer, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)
deer_v <- sf::st_as_sf(deer_v)

# Unique identifiers
un.id <- unique(deer$key)

sub_avail <- data.frame()

## Load grids for available points
# Remove available locations that are in the same cell as used points
for (i in 1:length(un.id)) {
  
  temp <- deer %>% filter(key==un.id[i])
  pts <- deer_v %>% filter(key==un.id[i])
  pts_u <- pts %>% filter(case==1)
  pts_a <- pts %>% filter(case==0)
  
  ## Load grids for available points
  filename <- paste0("output/deer_avail_cells/", un.id[i], "-avail.shp")
  grid <- vect(filename)
  
  # Extract the cell numbers that have a point in them
  grid_u <- extract(grid, pts_u)[,-1]
  
  # subset the polygon grid by only the cells that have points in them
  grid_a <- grid %>% filter(!(cell_no %in% grid_u))
  
  # Crop out cellsthat have used points in the 
  keep_a <- crop(pts_a, grid_a)
  
  # Extract RSF values
  a_extr <- extract(pred, keep_a)
  
  # Add columns for ID and case
  a_extr$key <- un.id[i]
  a_extr$case <- 0
  
  # Save to data frame
  sub_avail <- bind_rows(sub_avail, a_extr)
}

write_csv(sub_avail, "output/deer_avail_not_overlapping_used.csv")

# rearrange columns and convert to tibble
sub_avail <- sub_avail %>%
  select(key, case, RSF) %>%
  as_tibble()

# subset to used locations
used <- deer %>% filter(case==1)

# convert locations to spatvector
used <- vect(used, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)         

# Extract RSF values at each deer and available point
used_extr <- extract(pred, used)
used$RSF <- used_extr$RSF

# convert to data frame (will drop the X, Y coords, but thats ok)
used <- as.data.frame(used)
# drop extra columns and convert to tibble
used <- used %>% 
  select(key, case, RSF) %>% 
  as_tibble()

# How many points per key for used?
nobs <- used %>% group_by(key) %>% count()
un.id <- unique(used$key)

# 
avail <- data.frame()
for (i in 1:length(un.id)) {
  
  temp <- sub_avail %>% filter(key==un.id[i])
  n <- unname(unlist(as.vector(nobs[nobs$key==un.id[i],2])))
  
  temp <- temp %>% slice_sample(n=n)
  avail <- bind_rows(avail, temp)
}



## combine
dat <- bind_rows(used, avail)

# ggplot(data=deer) +
#   geom_histogram(aes(x=RSF)) +
#   facet_wrap(vars(case))
# 
quantile(dat$RSF[dat$case==1], probs=c(0.1,0.25, 0.5, 0.75, 1))

## Calculate precision and recall (sensitivity) for case==1
# Create vector of thresholds
thresh <- seq(0.01,1,by=0.01)
res <- matrix(NA, ncol=3, nrow=length(thresh))
res[,1] <- thresh
for (i in 1:length(thresh)) {
  
  # Convert presence/available to factor
  actu <- factor(dat$case, levels=c(0,1))
  # Convert RSF to binary based on threshold
  class <- factor(if_else(dat$RSF >= thresh[i], 1, 0), levels=c(0,1))
  
  cm <- caret::confusionMatrix(class, actu, positive="1")
  
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

res <- res %>%
  mutate(Fmeas = 2*((Precision * Recall)/(Precision + Recall)))


ggplot(res, aes(x=threshold, y=Fmeas)) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  geom_point() +
  theme_bw()

# Area unter the Prcision-Recall curve
pracma::trapz(sort(res$Recall), res$Precision)

#
pracma::trapz(res$threshold, res$Fmeas)




m <- c(0, 0.247, 0,
       0.247, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc1 <- classify(deer_pred, rclmat, include.lowest=TRUE)
plot(rc1)

