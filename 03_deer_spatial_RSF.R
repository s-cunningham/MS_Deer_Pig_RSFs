library(tidyverse)
library(terra)
library(lme4)

theme_set(theme_bw())

#### Read Data ####
# locations (used & available)
deer <- read_csv("output/deer_used_avail_locations.csv") %>%
  # Add column for weight
  mutate(weight=if_else(case==1, 1, 5000))

# Rasters
rast_list <- c("data/landscape_data/bottomlandHW_180m_sum.tif",
               "data/landscape_data/decidmixed_180m_sum.tif",
               "data/landscape_data/othercrops_180m_sum.tif",
               "data/landscape_data/gramanoids_180m_sum.tif", 
               "data/landscape_data/evergreen_180m_sum.tif",
               "data/landscape_data/barren_180m_sum.tif", 
               "data/landscape_data/developed_180m_sum.tif",
               "data/landscape_data/palatable_crops_180m_sum.tif",
               "data/landscape_data/water_180m_sum.tif") # replace with water (distance)
layers <- rast(rast_list)

# Reclassify missing data to 0
m <- rbind(c(NA, 0))
layers <- classify(layers, m)

# Convert to % 
layers <- layers / 149

# Center and scale continuous rasters
layers <- scale(layers)

# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(layers))

# convert locations to spatvector
deer_v <- vect(deer, geom=c("X", "Y"), crs=crs(layers), keepgeom=FALSE)

# Rename layers
names(layers) <- c("bottomlandhw", "decidmixed", "othercrops", "gramanoids", "evergreen", "barren", "developed", "foodcrops", "water")

#### Extract covariates and used and available locations ####

# Extract distance values at used and available locations
dat_deer <- extract(layers, deer_v)

# Join extracted data back to location data frame
deer <- bind_cols(deer, dat_deer)

ggplot(deer,aes(x=water, y=case)) +
  geom_point(alpha=0.3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) 

ggplot(deer, aes(x=bottomlandhw, y=foodcrops, color=id)) +
  geom_point(alpha=0.1) +
  facet_wrap(vars(case)) +
  theme(legend.position="none")


cor(deer[,c(8:16)])
