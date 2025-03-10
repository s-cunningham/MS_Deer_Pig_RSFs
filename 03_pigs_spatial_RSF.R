library(tidyverse)
library(terra)

#### Read Data ####
# locations (used & available)
pigs <- read_csv("output/pigs_used_avail_locations.csv")

# mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 

# Rasters
rast_list <- c("bottomlandHW_210m_sum.tif",
               "decidmixed_210m_sum.tif",
               "othercrops_210m_sum.tif",
               "gramanoids_210m_sum.tif", 
               "evergreen_210m_sum.tif",
               "barren_210m_sum.tif", 
               "developed_210m_sum.tif",
               "palatable_crops_210m_sum.tif") # water (distance)
layers <- rast(rast_list)


