library(tidyverse)
library(terra)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
deer <- read_csv("output/deer_used_avail_locations.csv") %>%
  # Add column for weight
  mutate(weight=if_else(case==1, 1, 5000))

deer %>% group_by(key, case) %>% count() %>% 
  pivot_wider(names_from="case", values_from="n") %>%
  mutate(ratioUA=`1`/`0`) %>% 
  ggplot() + geom_density(aes(x=ratioUA))

# Rasters
rast_list <- c("data/landscape_data/allhardwoods_180m_sum.tif",
               "data/landscape_data/shrublands_180m_sum.tif",
               "data/landscape_data/gramanoids_180m_sum.tif", 
               "data/landscape_data/palatable_crops_180m_sum.tif",
               "data/landscape_data/developed_180m_sum.tif") 

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
names(layers) <- c("allhardwoods", "shrubs", "gramanoids", "foodcrops", "developed", "water_dist")

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

# put the layer names back
names(layers) <- c("allhardwoods", "shrubs", "gramanoids", "foodcrops", "developed", "water_dist")

# crop to map extent
layers <- crop(layers, ms)

# Get mean and sd from entire covariate rasters
mean_vals <- global(layers, "mean", na.rm = TRUE)
sd_vals   <- global(layers, "sd", na.rm = TRUE)

# Apply to covariate rasters
cov_scaled <- (layers - mean_vals$mean) / sd_vals$sd

rast_cs <- bind_cols(mean_vals, sd_vals)
rast_cs <- rownames_to_column(rast_cs, var="layer")
write_csv(rast_cs, "output/deer_raster_mean_sds.csv")

#### Extract covariates and used and available locations ####

## Reformat shapefiles
# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(layers))

# convert locations to spatvector
deer_v <- vect(deer, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)

# Extract distance values at used and available locations
dat_deer <- extract(layers, deer_v)

# Join extracted data back to location data frame
deer <- bind_cols(deer, dat_deer)

# correlation matrix
cor(deer[,c(7:19)])

# write file so we don't always have to wait for the rasters to do stuff
write_csv(deer, "output/deer_used_avail_covariates.csv")
