# ----------------------------------------------------------------------------------------------------------------------------
# Name:
# Author: Stephanie Cunningham
# Objective: Build a geographically weighted raster for each of the four pig studies
# Input: Pig locations (csv), raster to extract CRS and extent
# Output:
# ----------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(terra)

aea_crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

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

# Read in raster as a template for CRS
forest <- rast("data/landscape_data/allhardwoodsMODIFIED_180m_sum.tif")

# Read in pigs location file
pigs <- read_csv("data/location_data/pigs_filtered.csv")

# Convert to SpatVector
pigs <- vect(pigs, geom=c("X", "Y"), crs=crs(forest))

# Find centroid of locations
pigs <- split(pigs, "study")

# Calculate mean point
pigs <- lapply(pigs, find_study_area_center)



# Extract extent from raster
ms_ext <- ext(forest)
ms_crs <- crs(forest, proj=TRUE)

r <- rast(ms_ext, resolution=30)
crs(r) <- crs(forest)

test <- distance(r, pigs$Delta)

## --- Building an geographically weighted raster --- ##

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
  
  # Normalize weights to sum to 1 if needed
  # weights <- weights / sum(weights)
  
  # Assign weights to the raster
  values(r) <- weights
  
  # Return raster
  return(r)
  
}


# Create rasters
gwr <- lapply(pigs, dist_weight_raster, extent=ms_ext)




