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
rast_list <- c("data/landscape_data/evergreen_180m_sum.tif",
               "data/landscape_data/deciduous_180m_sum.tif",
               "data/landscape_data/mixed_180m_sum.tif",
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
# water <- rast("data/landscape_data/RSinterarealMerge_distance30m.tif")
# water <- resample(water, layers)
# ext(water) <- ext(layers)

# layers <- c(layers, water)

# Center and scale continuous rasters
# layers <- scale(layers)

# Rename layers
names(layers) <- c("evergreen", "deciduous", "mixed", "shrubs", "othercrops",
                   "gramanoids", "bottomland", "herbwetl", "foodcrops", "developed", "water")

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
cor(deer[,c(7:17)])

# write file so we don't always have to wait for the rasters to do stuff
write_csv(deer, "output/deer_used_avail_covariates.csv")
