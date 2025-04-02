# ----------------------------------------------------------------------------------------------------------------------------
# Name:
# Author: Stephanie Cunningham
# Objective: Build a geographically weighted raster for each of the four pig studies
# Input: Pig locations (csv), raster to extract CRS and extent
# Output:
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(tidyterra)

# aea_crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Function to find the centroid of a group of points
find_study_area_center <- function(points_spatvector) {
  # Take just the coordinates of each point
  coords <- crds(points_spatvector)
  # Calculate the centroid (mean of coordinates)
  center_x <- mean(coords[,1])
  center_y <- mean(coords[,2])
  # Create a data frame
  df <- data.frame(X=center_x, Y=center_y)
  # Convert to SpatVect
  df <- vect(df, geom=c("X", "Y"), crs=crs(points_spatvector))
  
  return(df)
}

## Function to create a geographically weighted raster using terra
dist_weight_raster <- function(centroid, extent, coord_ref, resolution=30, distance_decay_factor=1) {
  # distance_decay_factor = 1: Linear inverse distance weighting
  # distance_decay_factor = 2: Quadratic inverse distance weighting (weights decrease more rapidly)
  
  # Create raster with specified extent and resolution (resolution defaults to 30)
  r <- rast(extent, resolution=resolution)
  
  # Define projection for empty raster
  crs(r) <- crs(coord_ref) 
  
  # Calculate distances from each raster cell to the centroid
  distances <- distance(r, centroid)
  
  ## Convert distances to inverse weights
  # Add small value to avoid division by 0
  weights <- 1 / (distances + 0.0001)^distance_decay_factor

  # mask and crop mississippi
  w <- crop(weights, ms_full, mask=TRUE)
  
  # Return raster
  return(w)
  
}


## --- Building a geographically weighted raster --- ##

# Read in raster as a template for CRS
forest <- rast("data/landscape_data/allhardwoodsMODIFIED_180m_sum.tif")

# Read in Mississippi shapefile
ms_full <- vect("data/landscape_data/mississippi_ACEA.shp") 
ms_full <- project(ms_full, crs(forest))

# Read in pigs location file
pigs <- read_csv("data/location_data/pigs_filtered.csv")

# Convert to SpatVector
# Note that these locations are in UTM coordinates, and need to be reprojected to AEA to match the working projection
pigs <- vect(pigs, geom=c("X", "Y"), crs="EPSG:32616")
pigs <- project(pigs, crs(forest))

# Find centroid of locations
pigs <- split(pigs, "study")

# Calculate mean point
pigs <- lapply(pigs, find_study_area_center)

# Extract extent from raster
ms_ext <- ext(forest)
ms_crs <- crs(forest, proj=TRUE)

# Create rasters
gwr <- lapply(pigs, dist_weight_raster, extent=ms_ext, coord_ref=ms_crs)


r <- rast(ms_ext, resolution=30)
# Define projection for empty raster
crs(r) <- crs(ms_crs) 

# Calculate distances from each raster cell to the centroid
distances <- distance(r, pigs$Delta)

## Convert distances to inverse weights
# Add small value to avoid division by 0
weights <- 1 / (distances + 0.0001)^distance_decay_factor

# Normalize weights to sum to 1 if needed
# weights <- weights / sum(weights)

# Assign weights to the raster
values(r) <- weights

ggplot() +
  geom_spatraster(data=forest) +
  geom_spatvector(data=pigs$Delta, color="red", size=3)

test <- pigs$Delta